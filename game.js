
const SKY = [0.62, 0.78, 0.92];
const SAND = [0.83, 0.76, 0.6];
const SEA = [0.16, 0.5, 0.66];
const ROADC = [[0.24,0.25,0.28],[0.25,0.26,0.29],[0.27,0.28,0.31],[0.3,0.31,0.34],[0.33,0.34,0.37],[0.36,0.36,0.39],[0.4,0.4,0.42]];

async function loadData() {
  const r = await fetch('data.gz'); if (!r.ok) throw new Error('map data missing');
  const ds = r.body.pipeThrough(new DecompressionStream('gzip'));
  return JSON.parse(await new Response(ds).text());
}

// ---------- geometry helpers ----------
function earcut(P) { // P: [[x,z]...] CCW in x,z plane; returns index triples
  const n = P.length, idx = [...Array(n).keys()], out = [];
  let guard = 0;
  const cross = (a, b, c) => (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);
  while (idx.length > 3 && guard++ < 4000) {
    let cut = false;
    for (let i = 0; i < idx.length; i++) {
      const a = P[idx[(i+idx.length-1)%idx.length]], b = P[idx[i]], c = P[idx[(i+1)%idx.length]];
      if (cross(a,b,c) <= 1e-9) continue;
      let ok = true;
      for (let j = 0; j < idx.length; j++) {
        const p = P[idx[j]];
        if (p === a || p === b || p === c) continue;
        if (cross(a,b,p) >= 0 && cross(b,c,p) >= 0 && cross(c,a,p) >= 0) { ok = false; break; }
      }
      if (ok) { out.push(idx[(i+idx.length-1)%idx.length], idx[i], idx[(i+1)%idx.length]); idx.splice(i,1); cut = true; break; }
    }
    if (!cut) break;
  }
  if (idx.length === 3) out.push(idx[0], idx[1], idx[2]);
  else if (idx.length > 3) for (let i = 1; i < idx.length-1; i++) out.push(idx[0], idx[i], idx[i+1]);
  return out;
}
class MB { // mesh builder
  constructor() { this.v = []; }
  tri(a, b, c, n, col) { for (const p of [a,b,c]) this.v.push(p[0],p[1],p[2],n[0],n[1],n[2],col[0],col[1],col[2]); }
  quad(a,b,c,d,n,col) { this.tri(a,b,c,n,col); this.tri(a,c,d,n,col); }
  norm(a,b,c){const ux=b[0]-a[0],uy=b[1]-a[1],uz=b[2]-a[2],vx=c[0]-a[0],vy=c[1]-a[1],vz=c[2]-a[2];let x=uy*vz-uz*vy,y=uz*vx-ux*vz,z=ux*vy-uy*vx;const l=Math.hypot(x,y,z)||1;return[x/l,y/l,z/l];}
  box(cx,cy,cz,sx,sy,sz,col,rot=0){ // centred in xz, base at cy
    const c=Math.cos(rot),s=Math.sin(rot);
    const P=(x,y,z)=>[cx+x*c-z*s,cy+y,cz+x*s+z*c];
    const hx=sx/2,hz=sz/2,h=sy;
    const A=P(-hx,0,-hz),B=P(hx,0,-hz),C=P(hx,0,hz),D=P(-hx,0,hz),E=P(-hx,h,-hz),F=P(hx,h,-hz),G=P(hx,h,hz),H=P(-hx,h,hz);
    this.quad(E,H,G,F,[0,1,0],col);
    this.quad(A,B,F,E,this.norm(A,B,F),shade(col,.8));
    this.quad(B,C,G,F,this.norm(B,C,G),shade(col,.9));
    this.quad(C,D,H,G,this.norm(C,D,H),shade(col,.8));
    this.quad(D,A,E,H,this.norm(D,A,E),shade(col,.9));
  }
  cyl(cx,cy,cz,r,h,col,seg=10,r2=null){
    r2=r2==null?r:r2;
    for(let i=0;i<seg;i++){const a0=i/seg*6.2832,a1=(i+1)/seg*6.2832;
      const A=[cx+Math.cos(a0)*r,cy,cz+Math.sin(a0)*r],B=[cx+Math.cos(a1)*r,cy,cz+Math.sin(a1)*r];
      const C=[cx+Math.cos(a1)*r2,cy+h,cz+Math.sin(a1)*r2],D=[cx+Math.cos(a0)*r2,cy+h,cz+Math.sin(a0)*r2];
      this.quad(A,B,C,D,this.norm(A,B,C),shade(col,0.75+0.25*Math.abs(Math.cos((a0+a1)/2))));
      if(r2<1e-3) {} else this.tri(D,C,[cx,cy+h,cz],[0,1,0],col);}
  }
  dome(cx,cy,cz,r,col,seg=14,rings=5){
    for(let k=0;k<rings;k++){const p0=k/rings*1.5708,p1=(k+1)/rings*1.5708;
      for(let i=0;i<seg;i++){const a0=i/seg*6.2832,a1=(i+1)/seg*6.2832;
        const f=(p,a)=>[cx+Math.cos(a)*Math.cos(p)*r,cy+Math.sin(p)*r,cz+Math.sin(a)*Math.cos(p)*r];
        const A=f(p0,a0),B=f(p0,a1),C=f(p1,a1),D=f(p1,a0);
        this.quad(A,B,C,D,this.norm(A,B,C),shade(col,0.8+0.2*Math.sin(p0+0.3)));}}
  }
}
function shade(c,f){return[c[0]*f,c[1]*f,c[2]*f];}
function rnd(i){const x=Math.sin(i*127.1+311.7)*43758.5453;return x-Math.floor(x);}

const BCOL=[[0.93,0.91,0.86],[0.9,0.87,0.8],[0.96,0.94,0.9],[0.88,0.84,0.76],[0.94,0.9,0.82],[0.86,0.86,0.86]];

async function startGame(root, setStatus) {
  setStatus('Unpacking Muscat map...');
  const D = await loadData();
  const T = D.T;
  const tiles = new Map();
  for (const t of D.t) tiles.set(t[0]+','+t[1], { tx:t[0], tz:t[1], r:t[2], b:t[3], m:t[4], c:t[5], mesh:null, polys:null });
  const pois = D.p;

  // DOM
  root.innerHTML = '';
  if (getComputedStyle(root).position === 'static') root.style.position = 'relative'; root.style.overflow = 'hidden'; root.style.background = '#a0c8ea'; root.style.touchAction = 'none'; root.style.userSelect='none';
  const cv = document.createElement('canvas'); cv.style.cssText = 'width:100%;height:100%;display:block'; root.appendChild(cv);
  const hud = document.createElement('div');
  hud.style.cssText = 'position:absolute;inset:0;pointer-events:none;font:600 14px system-ui,sans-serif;color:#fff;text-shadow:0 1px 3px #000';
  hud.innerHTML = `
   <div id="tx-top" style="position:absolute;left:10px;top:10px;display:flex;flex-direction:column;gap:6px">
     <div id="tx-money" style="font-size:20px">0.000 OMR</div>
     <div style="width:120px;height:10px;background:#0006;border-radius:5px;overflow:hidden"><div id="tx-fuel" style="height:100%;width:100%;background:#4ade80"></div></div>
     <div style="font-size:11px;opacity:.9">FUEL</div>
   </div>
   <div id="tx-msg" style="position:absolute;left:0;right:0;top:12px;text-align:center;font-size:15px;padding:0 130px"></div>
   <div id="tx-spd" style="position:absolute;right:12px;bottom:12px;font-size:28px;text-align:right">0<span style="font-size:12px"> km/h</span></div>
   <div id="tx-place" style="position:absolute;left:10px;bottom:10px;font-size:12px;opacity:.9"></div>
   <div id="tx-arrow" style="position:absolute;left:50%;top:76px;margin-left:-18px;width:36px;height:36px;font-size:30px;text-align:center;line-height:36px;transition:none">▲</div>
   <div id="tx-dist" style="position:absolute;left:50%;top:112px;transform:translateX(-50%);font-size:12px"></div>`;
  root.appendChild(hud);
  const mini = document.createElement('canvas'); mini.width = 150; mini.height = 150;
  mini.style.cssText = 'position:absolute;right:10px;top:10px;width:110px;height:110px;border-radius:50%;border:2px solid #fff9;background:#0007';
  root.appendChild(mini);
  const mctx = mini.getContext('2d');
  const msg = hud.querySelector('#tx-msg'), moneyEl = hud.querySelector('#tx-money'), fuelEl = hud.querySelector('#tx-fuel'),
    spdEl = hud.querySelector('#tx-spd'), arrow = hud.querySelector('#tx-arrow'), distEl = hud.querySelector('#tx-dist'), placeEl = hud.querySelector('#tx-place');

  // touch controls
  const keys = {};
  const btn = (txt, css, k) => { const b = document.createElement('div'); b.textContent = txt; b.style.cssText = 'position:absolute;width:64px;height:64px;border-radius:32px;background:#0008;color:#fff;font:700 26px system-ui;display:flex;align-items:center;justify-content:center;border:2px solid #fff6;'+css;
    const on = e => { keys[k] = true; e.preventDefault(); }, off = e => { keys[k] = false; e.preventDefault(); };
    b.addEventListener('pointerdown', on); b.addEventListener('pointerup', off); b.addEventListener('pointerleave', off); b.addEventListener('pointercancel', off); root.appendChild(b); return b; };
  const touch = matchMedia('(pointer:coarse)').matches || 'ontouchstart' in window;
  if (touch) {
    btn('◀','left:14px;bottom:96px','left'); btn('▶','left:92px;bottom:96px','right');
    btn('GO','right:14px;bottom:76px;width:76px;height:76px;border-radius:38px;background:#16a34a99','up'); btn('■','right:104px;bottom:96px;font-size:20px','brake');
  }
  addEventListener('keydown', e => { const m = {ArrowUp:'up',KeyW:'up',ArrowDown:'brake',KeyS:'brake',ArrowLeft:'left',KeyA:'left',ArrowRight:'right',KeyD:'right',Space:'brake',KeyR:'reset'}[e.code]; if (m) { keys[m] = true; e.preventDefault(); } });
  addEventListener('keyup', e => { const m = {ArrowUp:'up',KeyW:'up',ArrowDown:'brake',KeyS:'brake',ArrowLeft:'left',KeyA:'left',ArrowRight:'right',KeyD:'right',Space:'brake',KeyR:'reset'}[e.code]; if (m) keys[m] = false; });

  // GL
  const gl = cv.getContext('webgl2', { antialias: true });
  if (!gl) { setStatus('This browser has no WebGL2.'); return; }
  const vs = `#version 300 es
  in vec3 p; in vec3 n; in vec3 c; uniform mat4 vp; uniform vec3 cam; uniform mat4 mdl; out vec3 vc; out vec3 wq;
  void main(){ vec3 L=normalize(vec3(.4,.9,.3)); float d=abs(dot(n,L))*.5+.55; vc=c*d; vec4 wp=mdl*vec4(p,1.); gl_Position=vp*wp; wq=wp.xyz; }`;
  const fs = `#version 300 es
  precision highp float; in vec3 vc; in vec3 wq; uniform vec3 sky; uniform vec3 cam; out vec4 o;
  void main(){ float fg=length(wq-cam); float f=clamp((fg-350.)/950.,0.,1.); o=vec4(mix(vc,sky,f),1.); }`;
  const sh = (t, s) => { const o = gl.createShader(t); gl.shaderSource(o, s); gl.compileShader(o); if (!gl.getShaderParameter(o, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(o)); return o; };
  const prog = gl.createProgram(); gl.attachShader(prog, sh(gl.VERTEX_SHADER, vs)); gl.attachShader(prog, sh(gl.FRAGMENT_SHADER, fs)); gl.linkProgram(prog); gl.useProgram(prog);
  const aP = gl.getAttribLocation(prog,'p'), aN = gl.getAttribLocation(prog,'n'), aC = gl.getAttribLocation(prog,'c');
  const uMdl = gl.getUniformLocation(prog,'mdl'), uVP = gl.getUniformLocation(prog,'vp'), uCam = gl.getUniformLocation(prog,'cam'), uSky = gl.getUniformLocation(prog,'sky');
  gl.enable(gl.DEPTH_TEST); 
  function upload(arr) {
    const vao = gl.createVertexArray(); gl.bindVertexArray(vao);
    const buf = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, buf); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(arr), gl.STATIC_DRAW);
    for (const [a, o] of [[aP,0],[aN,3],[aC,6]]) { gl.enableVertexAttribArray(a); gl.vertexAttribPointer(a, 3, gl.FLOAT, false, 36, o*4); }
    return { vao, buf, count: arr.length / 9 };
  }
  const freeMesh = m => { if (m) { gl.deleteVertexArray(m.vao); gl.deleteBuffer(m.buf); } };

  // tile building
  function buildTile(t) {
    const mb = new MB(); const ox = t.tx*T, oz = t.tz*T;
    // sea cells (50 m raster runs: row, col0, col1)
    for (const c of t.c) { const z0 = oz + c[0]*50, x0 = ox + c[1]*50, x1 = ox + (c[2]+1)*50, z1 = z0 + 50;
      mb.quad([x0,0.02,z0],[x0,0.02,z1],[x1,0.02,z1],[x1,0.02,z0],[0,1,0],SEA); }
    // mountains
    for (let i = 0; i < t.m.length; i++) { const m = t.m[i]; const x=ox+m[0]/10, z=oz+m[1]/10, r=m[2], h=m[3];
      const col = [0.58+rnd(i)*0.1, 0.46+rnd(i+3)*0.08, 0.38+rnd(i+5)*0.06]; mb.cyl(x,0,z,r,h,col,7,r*0.12); }
    // roads
    for (const r of t.r) {
      const w = r[1], col = ROADC[r[0]], p = r[3], y = 0.05 + (6-r[0])*0.004;
      for (let i = 0; i+3 < p.length; i += 2) {
        const x0 = ox+p[i]/10, z0 = oz+p[i+1]/10, x1 = ox+p[i+2]/10, z1 = oz+p[i+3]/10;
        const dx = x1-x0, dz = z1-z0, l = Math.hypot(dx,dz); if (l < 0.01) continue;
        const nx = -dz/l*w/2, nz = dx/l*w/2;
        const e = 0.6, ex = dx/l*e, ez = dz/l*e;
        const A=[x0-ex+nx,y,z0-ez+nz],B=[x1+ex+nx,y,z1+ez+nz],C=[x1+ex-nx,y,z1+ez-nz],Dd=[x0-ex-nx,y,z0-ez-nz];
        mb.quad(A,Dd,C,B,[0,1,0],col); mb.quad(A,B,C,Dd,[0,1,0],col);
        if (r[0] <= 3 && w >= 9 && ((i/2)%2===0)) { // lane dashes
          const q=0.3,cx0=x0+dx*0.3,cz0=z0+dz*0.3,cx1=x0+dx*0.6,cz1=z0+dz*0.6;
          const a=[cx0+nx*.03,y+.012,cz0+nz*.03],b=[cx1+nx*.03,y+.012,cz1+nz*.03],c=[cx1-nx*.03,y+.012,cz1-nz*.03],d=[cx0-nx*.03,y+.012,cz0-nz*.03];
          mb.quad(a,d,c,b,[0,1,0],[0.9,0.9,0.85]); mb.quad(a,b,c,d,[0,1,0],[0.9,0.9,0.85]);
        }
      }
    }
    // buildings
    const polys = [];
    t.b.forEach((b, bi) => {
      const h = b[0]/10, p = b[1]; const n = p.length/2; let P = [];
      for (let i = 0; i < n; i++) P.push([ox+p[2*i]/10, oz+p[2*i+1]/10]);
      let a2 = 0; for (let i = 0; i < n; i++) { const q = P[(i+1)%n]; a2 += P[i][0]*q[1]-q[0]*P[i][1]; }
      if (a2 > 0) P.reverse(); // want CCW seen from above in x,z with z south... normalise so wall winding below is consistent
      const col = BCOL[(bi*7+t.tx*3+t.tz)%BCOL.length];
      const roofc = shade(col, 0.82);
      let minx=1e9,maxx=-1e9,minz=1e9,maxz=-1e9; for (const q of P) { minx=Math.min(minx,q[0]);maxx=Math.max(maxx,q[0]);minz=Math.min(minz,q[1]);maxz=Math.max(maxz,q[1]); }
      polys.push({P, minx, maxx, minz, maxz});
      for (let i = 0; i < n; i++) {
        const q0 = P[i], q1 = P[(i+1)%n];
        const A=[q0[0],0,q0[1]],B=[q1[0],0,q1[1]],C=[q1[0],h,q1[1]],Dd=[q0[0],h,q0[1]];
        const nn = mb.norm(A,B,C); const f = 0.78+0.22*Math.abs(nn[0]);
        mb.quad(A,B,C,Dd,nn,shade(col,f)); mb.quad(A,Dd,C,B,[-nn[0],-nn[1],-nn[2]],shade(col,f));
        const el = Math.hypot(q1[0]-q0[0], q1[1]-q0[1]);
        if (h > 5 && el > 5) { const ux = (q1[0]-q0[0])/el, uz = (q1[1]-q0[1])/el, ox2 = nn[0]*0.06, oz2 = nn[2]*0.06, wc = [0.22,0.3,0.38];
          const floors = Math.min(14, Math.floor((h-1)/3.2));
          for (let fl = 0; fl < floors; fl++) { const y0 = 1.3 + fl*3.2, y1 = y0 + 1.3;
            const a0=[q0[0]+ux*1.2+ox2,y0,q0[1]+uz*1.2+oz2], b0=[q1[0]-ux*1.2+ox2,y0,q1[1]-uz*1.2+oz2], c0=[b0[0],y1,b0[2]], d0=[a0[0],y1,a0[2]];
            mb.quad(a0,b0,c0,d0,nn,wc); const a1=[a0[0]-2*ox2,y0,a0[2]-2*oz2], b1=[b0[0]-2*ox2,y0,b0[2]-2*oz2]; mb.quad(a1,[a1[0],y1,a1[2]],[b1[0],y1,b1[2]],b1,[-nn[0],-nn[1],-nn[2]],wc); } }
      }
      const idx = earcut(P);
      for (let i = 0; i < idx.length; i += 3) { const a=P[idx[i]],b=P[idx[i+1]],c=P[idx[i+2]];
        const A=[a[0],h,a[1]],B=[b[0],h,b[1]],C=[c[0],h,c[1]]; mb.tri(A,B,C,[0,1,0],roofc); mb.tri(A,C,B,[0,1,0],roofc); }
    });
    // landmarks
    for (const p of pois) {
      if (Math.floor(p[0]/T) !== t.tx || Math.floor(p[1]/T) !== t.tz || p[2] !== 'p') continue;
      const nm = p[3];
      if (/Grand Mosque/.test(nm)) { mb.box(p[0],0,p[1],60,18,60,[0.97,0.95,0.9]); mb.dome(p[0],18,p[1],22,[0.9,0.8,0.45],16,6);
        for (const [dx,dz] of [[-48,-48],[48,-48],[-48,48],[48,48]]) { mb.cyl(p[0]+dx,0,p[1]+dz,2.6,50,[0.97,0.95,0.9],10); mb.cyl(p[0]+dx,50,p[1]+dz,3.4,2,[0.9,0.8,0.45],10); mb.cyl(p[0]+dx,52,p[1]+dz,2,6,[0.97,0.95,0.9],10,0.1); } }
      else if (/Fort|Mirani/.test(nm)) { mb.cyl(p[0],0,p[1],11,22,[0.78,0.64,0.45],14); mb.cyl(p[0],22,p[1],12,3,[0.7,0.56,0.4],14); }
      else if (/Opera/.test(nm)) { mb.cyl(p[0]+0,0,p[1],5,40,[0.95,0.92,0.85],8); mb.dome(p[0],40,p[1],8,[0.8,0.7,0.45],10,4); }
    }
    t.polys = polys;
    t.mesh = upload(mb.v);
  }
  function ensure(cx, cz) {
    const R = 2; const keep = new Set();
    for (let i = -R; i <= R; i++) for (let j = -R; j <= R; j++) {
      const k = (cx+i)+','+(cz+j); keep.add(k); const t = tiles.get(k); if (t && !t.mesh) { buildTile(t); return; } // one per frame
    }
    for (const [k, t] of tiles) if (t.mesh && !keep.has(k)) { const dx=Math.abs(t.tx-cx),dz=Math.abs(t.tz-cz); if (dx>R+1||dz>R+1) { freeMesh(t.mesh); t.mesh=null; t.polys=null; } }
  }
  // taxi mesh
  const tb = new MB(); const Y=[0.98,0.78,0.08], DK=[0.12,0.12,0.14];
  tb.box(0,0.35,0,1.9,0.7,4.3,Y); tb.box(0,1.05,-0.2,1.7,0.6,2.2,[0.95,0.97,1]); tb.box(0,1.62,-0.2,1.78,0.1,2.3,Y); tb.box(0,1.72,-0.2,0.7,0.22,0.3,[1,1,1]);
  tb.box(-0.95,0,1.35,0.3,0.55,0.8,DK); tb.box(0.95,0,1.35,0.3,0.55,0.8,DK); tb.box(-0.95,0,-1.35,0.3,0.55,0.8,DK); tb.box(0.95,0,-1.35,0.3,0.55,0.8,DK);
  tb.box(0,0.4,2.15,1.6,0.14,0.05,[1,1,0.7]);
  const taxiMesh = upload(tb.v);
  // markers
  const mk = (col) => { const m = new MB(); m.cyl(0,0,0,2.2,60,col,10,2.2); return upload(m.v); };
  const pickM = mk([0.2,0.9,0.4]), dropM = mk([0.2,0.55,1]), fuelM = mk([1,0.45,0.1]);

  // math
  const persp = (fov, asp, n, f) => { const t = 1/Math.tan(fov/2); return [t/asp,0,0,0, 0,t,0,0, 0,0,(f+n)/(n-f),-1, 0,0,2*f*n/(n-f),0]; };
  const look = (e, c) => { let zx=e[0]-c[0], zy=e[1]-c[1], zz=e[2]-c[2]; let l=Math.hypot(zx,zy,zz); zx/=l;zy/=l;zz/=l;
    let xx=zz, xy=0, xz=-zx; l=Math.hypot(xx,xz); xx/=l; xz/=l; const yx=zy*xz-zz*xy, yy=zz*xx-zx*xz, yz=zx*xy-zy*xx;
    return [xx,yx,zx,0, xy,yy,zy,0, xz,yz,zz,0, -(xx*e[0]+xy*e[1]+xz*e[2]), -(yx*e[0]+yy*e[1]+yz*e[2]), -(zx*e[0]+zy*e[1]+zz*e[2]), 1]; };
  const mul = (a, b) => { const o = new Array(16).fill(0); for (let i=0;i<4;i++) for (let j=0;j<4;j++) for (let k=0;k<4;k++) o[j*4+i] += a[k*4+i]*b[j*4+k]; return o; };
  const trans = (x,y,z,rot) => { const c=Math.cos(rot),s=Math.sin(rot); return [c,0,-s,0, 0,1,0,0, s,0,c,0, x,y,z,1]; };
  // forward vector for heading h: (sin h, cos h) on x,z
  // road helpers
  const roadVerts = [];
  for (const t of D.t) for (const r of t[2]) if (r[0] <= 4 && r[0] >= 1) for (let i = 0; i+1 < r[3].length; i += 2) roadVerts.push([t[0]*T+r[3][i]/10, t[1]*T+r[3][i+1]/10, r]);
  const roadsAll = roadVerts;
  const rgrid = new Map();
  for (const t of D.t) for (const r of t[2]) { const p = r[3], hw = r[1]/2 + 3; for (let i = 0; i+3 < p.length; i += 2) {
    const x0 = t[0]*T+p[i]/10, z0 = t[1]*T+p[i+1]/10, x1 = t[0]*T+p[i+2]/10, z1 = t[1]*T+p[i+3]/10, l = Math.hypot(x1-x0, z1-z0), n = Math.max(1, Math.ceil(l/5));
    for (let k = 0; k <= n; k++) { const x = x0+(x1-x0)*k/n, z = z0+(z1-z0)*k/n, key = Math.floor(x/12)+','+Math.floor(z/12); let a = rgrid.get(key); if (!a) rgrid.set(key, a = []); a.push(x, z, hw); } } }
  const onRoad = (x, z) => { const gx = Math.floor(x/12), gz = Math.floor(z/12); for (let i = -1; i <= 1; i++) for (let j = -1; j <= 1; j++) { const a = rgrid.get((gx+i)+','+(gz+j)); if (a) for (let k = 0; k < a.length; k += 3) { const dx = a[k]-x, dz = a[k+1]-z; if (dx*dx+dz*dz < a[k+2]*a[k+2]) return true; } } return false; };
  const pickRoad = (cx, cz, dmin, dmax) => { const c = []; for (const v of roadsAll) { const d = Math.hypot(v[0]-cx, v[1]-cz); if (d >= dmin && d <= dmax) c.push(v); }
    if (c.length) return c[(Math.random()*c.length)|0]; let b = roadsAll[0], bd = 1e12; for (const v of roadsAll) { const d = Math.abs(Math.hypot(v[0]-cx, v[1]-cz) - dmin); if (d < bd) { bd = d; b = v; } } return b; };
  let start = pois.find(p => /Mutrah Fort/.test(p[3])) || pois[0];
  let sv = pickRoad(start[0], start[1], 0, 300);
  let best = 1e9; for (const v of roadsAll) { const d = Math.hypot(v[0]-start[0], v[1]-start[1]); if (d < best) { best = d; sv = v; } }
  const car = { x: sv[0], z: sv[1], h: 0, v: 0, fuel: 100, money: 0 };
  const startPos = { x: car.x, z: car.z };
  // heading along road
  { let bd=1e9; for (const t of D.t) for (const r of t[2]) { const p=r[3]; for (let i=0;i+3<p.length;i+=2){ const x=t[0]*T+p[i]/10,z=t[1]*T+p[i+1]/10; const d=Math.hypot(x-car.x,z-car.z); if(d<bd){bd=d;car.h=Math.atan2(t[0]*T+p[i+2]/10-x,t[1]*T+p[i+3]/10-z);} } } }
  const fuels = pois.filter(p => p[2] === 'f');
  window.__taxi = { car, get mission() { return mission; } };
  let mission = { state: 'wait', pick: null, drop: null, timer: 0, dist: 0, fare: 0 };
  let toast = '', toastT = 0;
  const say = (s, t=3) => { toast = s; toastT = t; };
  function newPickup() { const v = pickRoad(car.x, car.z, 250, 900); mission = { state: 'pick', pick: v, drop: null, timer: 0 }; }
  newPickup();

  function collides(x, z) {
    const cx = Math.floor(x/T), cz = Math.floor(z/T);
    for (let i = -1; i <= 1; i++) for (let j = -1; j <= 1; j++) { const t = tiles.get((cx+i)+','+(cz+j)); if (!t || !t.polys) continue;
      for (const b of t.polys) { if (x < b.minx || x > b.maxx || z < b.minz || z > b.maxz) continue;
        const P = b.P; let inside = false; for (let a = 0, c = P.length-1; a < P.length; c = a++) if (((P[a][1] > z) !== (P[c][1] > z)) && (x < (P[c][0]-P[a][0])*(z-P[a][1])/(P[c][1]-P[a][1])+P[a][0])) inside = !inside;
        if (inside) return true; } }
    return false;
  }
  function pointsOf(x, z, h) { const s=Math.sin(h), c=Math.cos(h), out=[]; for (const [a,b] of [[0.95,2.1],[-0.95,2.1],[0.95,-2.1],[-0.95,-2.1],[0,2.2]]) out.push([x + s*b + c*a, z + c*b - s*a]); return out; }

  // loop
  let last = performance.now(), camX = car.x, camZ = car.z, camH = car.h, alive = true, frame = 0, mmCache = null;
  const kmh = () => Math.abs(car.v)*3.6;
  function update(dt) {
    if (keys.reset) { car.x = startPos.x; car.z = startPos.z; car.v = 0; car.money = Math.max(0, car.money-2); car.fuel = Math.max(car.fuel, 30); say('Towed back to the start (-2 OMR)'); keys.reset = false; }
    const out = car.fuel <= 0;
    const acc = (keys.up && !out ? 11 : 0) - (keys.brake ? (car.v > 0.5 ? 26 : 8) : 0);
    car.v += acc*dt; car.v *= Math.pow(0.985, dt*60) ; if (!keys.up && !keys.brake) car.v *= Math.pow(0.99, dt*60);
    const road = onRoad(car.x, car.z); const maxv = road ? 34 : 9; if (!road && car.v > maxv) car.v *= Math.pow(0.9, dt*60); car.v = Math.max(-8, Math.min(maxv, car.v));
    const steer = (keys.left ? 1 : 0) - (keys.right ? 1 : 0);
    const turnRate = 1.7 * steer * (car.v / (Math.abs(car.v) + 4)) * (1 - Math.min(0.5, Math.abs(car.v)/maxv*0.5));
    car.h += turnRate*dt;
    const nx = car.x + Math.sin(car.h)*car.v*dt, nz = car.z + Math.cos(car.h)*car.v*dt;
    let hit = false; for (const p of pointsOf(nx, nz, car.h)) if (collides(p[0], p[1])) { hit = true; break; }
    if (hit) { car.v *= -0.25; } else { const d = Math.hypot(nx-car.x, nz-car.z); car.x = nx; car.z = nz; car.fuel -= d/450 + dt*0.02; }
    car.fuel = Math.max(0, car.fuel);
    if (out && Math.abs(car.v) < 0.5) say('Out of fuel. Press R (or tap ⛽ fuel station) - towing costs 2 OMR', 1);
    // missions
    const m = mission;
    if (m.state === 'pick') { const d = Math.hypot(m.pick[0]-car.x, m.pick[1]-car.z); m.dist = d;
      if (d < 9 && kmh() < 18) { const drop = pickRoad(car.x, car.z, 700, 2800); const dd = Math.hypot(drop[0]-car.x, drop[1]-car.z);
        mission = { state: 'drive', pick: null, drop, timer: dd/9 + 45, dist: dd, fare: 0.6 + dd/1000*0.45, d0: dd }; say('Passenger in. Take them to the blue beam.'); } }
    else if (m.state === 'drive') { m.timer -= dt; m.dist = Math.hypot(m.drop[0]-car.x, m.drop[1]-car.z);
      if (m.dist < 11 && kmh() < 20) { const bonus = Math.max(0, m.timer/60)*0.4; const pay = m.fare + bonus; car.money += pay; say('Fare ' + pay.toFixed(3) + ' OMR (tip ' + bonus.toFixed(3) + ')', 4); newPickup(); }
      else if (m.timer < -30) { say('Passenger got out. No fare.', 3); newPickup(); } }
    // fuel station
    for (const f of fuels) if (Math.hypot(f[0]-car.x, f[1]-car.z) < 14 && kmh() < 12 && car.fuel < 95) { const need = 100-car.fuel, cost = need*0.03; if (car.money >= cost) { car.money -= cost; car.fuel = 100; say('Refuelled (-' + cost.toFixed(3) + ' OMR)'); } else if (toastT < 0.5) say('Not enough money to refuel (' + cost.toFixed(3) + ' OMR)', 2); }
  }
  function drawMini() {
    const W = 150, c = mctx; c.clearRect(0,0,W,W); c.save(); c.translate(W/2, W/2); c.rotate(car.h + Math.PI); // heading up
    const sc = W/2/600; c.lineCap = 'round';
    const cx = Math.floor(car.x/T), cz = Math.floor(car.z/T);
    c.strokeStyle = '#ddd'; c.lineWidth = 1.4;
    for (let i=-1;i<=1;i++) for (let j=-1;j<=1;j++) { const t = tiles.get((cx+i)+','+(cz+j)); if (!t) continue; for (const r of t.r) { if (r[0] > 4) continue; const p = r[3]; c.beginPath(); for (let k=0;k+1<p.length;k+=2){ const x=(t.tx*T+p[k]/10-car.x)*sc, y=(t.tz*T+p[k+1]/10-car.z)*sc; k? c.lineTo(-x, y) : c.moveTo(-x, y);} c.stroke(); } }
    const dot = (v, col) => { let x = -(v[0]-car.x)*sc, y = (v[1]-car.z)*sc; const l = Math.hypot(x,y), R = 66; if (l > R) { x *= R/l; y *= R/l; } c.fillStyle = col; c.beginPath(); c.arc(x, y, 5, 0, 7); c.fill(); };
    if (mission.state === 'pick') dot(mission.pick, '#34d399'); if (mission.state === 'drive') dot(mission.drop, '#3b82f6');
    c.restore(); c.fillStyle = '#facc15'; c.beginPath(); c.moveTo(W/2, W/2-9); c.lineTo(W/2-6, W/2+7); c.lineTo(W/2+6, W/2+7); c.fill();
  }
  function draw(dt) {
    const w = cv.clientWidth, h = cv.clientHeight, dpr = Math.min(devicePixelRatio || 1, 2);
    if (cv.width !== (w*dpr|0) || cv.height !== (h*dpr|0)) { cv.width = w*dpr; cv.height = h*dpr; }
    gl.viewport(0,0,cv.width,cv.height); gl.clearColor(SKY[0],SKY[1],SKY[2],1); gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);
    const k = 1 - Math.pow(0.001, dt); camH += Math.atan2(Math.sin(car.h-camH), Math.cos(car.h-camH))*Math.min(1, dt*4); camX += (car.x-camX)*k; camZ += (car.z-camZ)*k;
    const back = (cv.width < cv.height ? 14 : 9) + kmh()*0.04, e = [car.x - Math.sin(camH)*back, 4.3 + kmh()*0.012, car.z - Math.cos(camH)*back], c = [car.x + Math.sin(camH)*6, 1.4, car.z + Math.cos(camH)*6];
    const vp = mul(persp(cv.width < cv.height ? 1.25 : 1.0, cv.width/cv.height, 0.5, 1600), look(e, c));
    const I4=new Float32Array([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]); gl.uniformMatrix4fv(uVP, false, new Float32Array(vp)); gl.uniformMatrix4fv(uMdl, false, I4); gl.uniform3fv(uCam, e); gl.uniform3fv(uSky, SKY);
    // ground
    const g = new MB(); const gx = Math.round(car.x/50)*50, gz = Math.round(car.z/50)*50, G = 1800;
    g.quad([gx-G,0,gz-G],[gx-G,0,gz+G],[gx+G,0,gz+G],[gx+G,0,gz-G],[0,1,0],SAND);
    const gm = upload(g.v); gl.bindVertexArray(gm.vao); gl.drawArrays(gl.TRIANGLES, 0, gm.count); freeMesh(gm);
    for (const t of tiles.values()) if (t.mesh) { const dx = (t.tx+0.5)*T - car.x, dz = (t.tz+0.5)*T - car.z; if (dx*dx+dz*dz > 2300*2300) continue; gl.bindVertexArray(t.mesh.vao); gl.drawArrays(gl.TRIANGLES, 0, t.mesh.count); }
    gl.uniformMatrix4fv(uMdl, false, new Float32Array(trans(car.x, 0, car.z, car.h)));
    gl.bindVertexArray(taxiMesh.vao); gl.drawArrays(gl.TRIANGLES, 0, taxiMesh.count);
    const marker = (v, m) => { gl.uniformMatrix4fv(uMdl, false, new Float32Array(trans(v[0], 0, v[1], 0))); gl.bindVertexArray(m.vao); gl.drawArrays(gl.TRIANGLES, 0, m.count); };
    if (mission.state === 'pick') marker(mission.pick, pickM); if (mission.state === 'drive') marker(mission.drop, dropM);
    for (const f of fuels) if (Math.hypot(f[0]-car.x, f[1]-car.z) < 1200) marker(f, fuelM);
  }
  function hudUpdate() {
    moneyEl.textContent = car.money.toFixed(3) + ' OMR'; fuelEl.style.width = car.fuel + '%'; fuelEl.style.background = car.fuel < 20 ? '#ef4444' : '#4ade80';
    spdEl.firstChild.textContent = Math.round(kmh());
    const m = mission; let t = '';
    if (toastT > 0) t = toast; else if (m.state === 'pick') t = 'Pick up the passenger (green beam)'; else if (m.state === 'drive') t = 'Drop-off ' + Math.round(m.dist) + ' m, ' + Math.max(0, Math.round(m.timer)) + 's for tip, fare ' + m.fare.toFixed(3) + ' OMR';
    msg.textContent = t;
    const tgt = m.state === 'pick' ? m.pick : m.state === 'drive' ? m.drop : null;
    if (tgt) { const a = Math.atan2(tgt[0]-car.x, tgt[1]-car.z); const rel = Math.atan2(Math.sin(a-car.h), Math.cos(a-car.h)); arrow.style.transform = 'rotate(' + (-rel) + 'rad)'; arrow.style.color = m.state === 'pick' ? '#34d399' : '#60a5fa'; distEl.textContent = Math.round(m.dist) + ' m'; }
    let near = ''; for (const p of pois) if (p[2] === 'p' && Math.hypot(p[0]-car.x, p[1]-car.z) < 400) { near = p[3]; break; } placeEl.textContent = near ? 'Near: ' + near : '';
  }
  setStatus('');
  function loop(now) {
    if (!alive) return; const dt = Math.min(0.05, (now-last)/1000); last = now; frame++;
    ensure(Math.floor(car.x/T), Math.floor(car.z/T)); update(dt); draw(dt); if (toastT > 0) toastT -= dt;
    if (frame % 3 === 0) { hudUpdate(); drawMini(); }
    requestAnimationFrame(loop);
  }
  // prebuild the starting tiles
  for (let i = 0; i < 12; i++) ensure(Math.floor(car.x/T), Math.floor(car.z/T));
  requestAnimationFrame(loop);
  return () => { alive = false; };
}
