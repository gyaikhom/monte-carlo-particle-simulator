\input boxes
vardef arrowhead expr p =
save A,u,a,b; pair A,u; path a,b;
A := point length(p) of p;
u := unitvector(direction length(p) of p);
a := A{-u}..(A - ahlength*u rotated 20);
b := A{-u}..(A - ahlength*u rotated -20);
( a & reverse(a) & b & reverse(b) ) --cycle
enddef;

outputtemplate := "neighbour-direct.mps";
beginfig(1);
  numeric u;
  u = 0.5cm;
  label.bot(btex Direct Mapping etex, (0.5[14,0], 0)*u);
  label.lft("1", (0, 1.5)*u);
  label.lft("0", (0, 2.5)*u);
  label.rt(btex $s$-field etex, (14.5, 4.5)*u);
  label.rt(btex $neighbour\_table$ etex, (14.5, 2.5)*u);
  fill (0,4)*u--(1,4)*u--(1,5)*u--(0,5)*u--cycle withcolor .85white;
  fill (3,4)*u--(4,4)*u--(4,5)*u--(3,5)*u--cycle withcolor .85white;
  fill (7,4)*u--(8,4)*u--(8,5)*u--(7,5)*u--cycle withcolor .85white;
  draw (0,0)*u--(0,3)*u--(14,3)*u--(14,0)*u;
  for r=1 upto 3: draw (0,r)*u--(14,r)*u; endfor
  draw (0,4)*u--(14,4)*u--(14,5)*u--(0,5)*u--cycle;
  for c=1 upto 8:
    draw (c,4)*u--(c,5)*u;
    label(decimal(c-1),(c-0.5,4.5)*u);
  endfor
  label(btex $\cdots$ etex,(9.5,4.5)*u);
  for c=11 upto 13:
    draw (c,4)*u--(c,5)*u;
    label(decimal(c + 29),(c+0.5,4.5)*u);
  endfor
  for c=1.5 upto 2.5:
    drawarrow (c,4)*u--(c,3)*u;
  endfor
  for c=4.5 upto 7:
    drawarrow (c,4)*u--(c,3)*u;
  endfor
  for c=11.5 upto 13.5:
    drawarrow (c,4)*u--(c,3)*u;
  endfor
endfig;

outputtemplate := "neighbour-tiered.mps";
beginfig(1);
  numeric u;
  u = 0.5cm;
  label.bot(btex Tiered Mapping etex, (0.5[11,0], 0)*u);
  label.lft("1", (0, 1.5)*u);
  label.lft("0", (0, 2.5)*u);
  label.rt(btex $s$-field etex, (14.5, 6.5)*u);
  label.rt(btex $neighbour\_idx\_table$ etex, (11.5, 4.5)*u);
  label.rt(btex $neighbour\_table$ etex, (11.5, 2.5)*u);
  draw (0,4)*u--(11,4)*u--(11,5)*u--(0,5)*u--cycle;
  for c=1 upto 5:
    draw (c,5)*u--(c,4)*u;
    drawarrow (c-0.5,4)*u--(c-0.5,3)*u;
    label(decimal(c-1),(c-0.5,4.5)*u);
  endfor
  label(btex $\cdots$ etex,(6.5,4.5)*u);
  for c=8 upto 10:
    draw (c,5)*u--(c,4)*u;
    drawarrow (c+0.5,4)*u--(c+0.5,3)*u;
    label(decimal(c+15),(c+0.5,4.5)*u);
  endfor
  draw (0,0)*u--(0,3)*u--(11,3)*u--(11,0)*u;
  for r=1 upto 3: draw (0,r)*u--(11,r)*u; endfor

  fill (0,6)*u--(1,6)*u--(1,7)*u--(0,7)*u--cycle withcolor .85white;
  fill (3,6)*u--(4,6)*u--(4,7)*u--(3,7)*u--cycle withcolor .85white;
  fill (7,6)*u--(8,6)*u--(8,7)*u--(7,7)*u--cycle withcolor .85white;
  draw (0,6)*u--(14,6)*u--(14,7)*u--(0,7)*u--cycle;
  for c=1 upto 8:
    draw (c,6)*u--(c,7)*u;
    label(decimal(c-1),(c-0.5,6.5)*u);
  endfor
  label(btex $\cdots$ etex,(9.5,6.5)*u);
  for c=11 upto 13:
    draw (c,6)*u--(c,7)*u;
    label(decimal(c + 29),(c+0.5,6.5)*u);
  endfor
  for c=1.5 upto 2.5:
    drawarrow (c,6)*u--(c-1,5)*u;
  endfor
  for c=4.5 upto 7:
    drawarrow (c,6)*u--(c-2,5)*u;
  endfor
  for c=11.5 upto 13.5:
    drawarrow (c,6)*u--(c-3,5)*u;
  endfor
endfig;

outputtemplate := "particle-repository.mps";
beginfig(1);
  numeric u;
  u = 0.5cm;
  for i=2 step 3 until 10:
    draw (i,5)*u--(i+2,5)*u--(i+2,1)*u--(i,1)*u--cycle;
    drawarrow (i-2,4.5)*u--(i,4.5)*u;
    for j=1 upto 4: draw (i,j)*u--(i+2,j)*u; endfor
  endfor
  draw (9,4.5)*u--(11,4.5)*u--(11,3.5)*u;
  draw (10.5,3.5)*u--(11.5,3.5)*u;
  draw (10.75,3.25)*u--(11.25,3.25)*u;
  draw (11,3)*u--(11.10,3)*u;
  label("8", (3, 3.5)*u);
  label("7", (3, 2.5)*u);
  label("6", (3, 1.5)*u);
  label("4", (6, 3.5)*u);
  label("3", (6, 2.5)*u);
  label("2", (6, 1.5)*u);
  label("5", (9, 3.5)*u);
  label("1", (9, 2.5)*u);

  picture tree;
  tree := image(
  path c[];
  pair p[];
  p[1] = (0,0);
  for i=2 upto 8:
    p[i] = p[i-1]+(u,0);
  endfor
  p[2] := p[2] shifted (0,1*u);
  p[3] := p[3] shifted (0,2*u);
  p[4] := p[4] shifted (0,1*u);
  p[5] := p[5] shifted (0,3*u);
  p[6] := p[6] shifted (0,1*u);
  p[7] := p[7] shifted (0,2*u);
  p[8] := p[8] shifted (0,1*u);
  c[1] = fullcircle scaled 1u;
  for i=2 upto 8:
    c[i] = c[1] shifted p[i];
    draw c[i];
  endfor
  draw c[1];
  label("1", p[1]);
  label("2", p[6]);
  label("3", p[4]);
  label("4", p[2]);
  label("5", p[8]);
  label("6", p[7]);
  label("7", p[3]);
  label("8", p[5]);
  draw p[5]--p[3] cutbefore c[5] cutafter c[3];
  draw p[5]--p[7] cutbefore c[5] cutafter c[7];
  draw p[3]--p[2] cutbefore c[3] cutafter c[2];
  draw p[3]--p[4] cutbefore c[3] cutafter c[4];
  draw p[7]--p[6] cutbefore c[7] cutafter c[6];
  draw p[7]--p[8] cutbefore c[7] cutafter c[8];
  draw p[2]--p[1] cutbefore c[2] cutafter c[1];
  );
  draw tree shifted (14u, 1.5u);
  pickup pencircle scaled 4pt;
  for i=0 step 3 until 9:
    drawdot (i,4.5)*u;
  endfor
endfig;

outputtemplate := "affine-merge.mps";
beginfig(1);
  numeric q;
  q := 0.75cm;
  circleit.u1(btex $\cup$ etex);
  circleit.u2(btex $\cup$ etex);
  circleit.i1(btex $\cap$ etex);
  circleit.i2(btex $\cap$ etex);
  circleit.i3(btex $\cap$ etex);
  circleit.t1(btex T etex);
  circleit.t2(btex T etex);
  circleit.r1(btex R etex);
  circleit.r2(btex R etex);
  circleit.s1(btex S etex);
  circleit.p0(btex 0 etex);
  circleit.p1(btex 1 etex);
  circleit.p2(btex 2 etex);
  circleit.p3(btex 3 etex);
  circleit.p4(btex 4 etex);
  circleit.p5(btex 5 etex);
  
  u1.c = (6*q, 4.75*q);
  u2.c = (5*q, 2.5*q);
  i1.c = (2*q, 1.75*q);
  i2.c = (4*q, 1.75*q);
  i3.c = (6*q, 1.75*q);
  t1.c = (5*q, 4*q);
  t2.c = (6*q, 3.25*q);
  r1.c = (7*q, 4*q);
  r2.c = (4*q, 3.25*q);
  s1.c = (3*q, 2.5*q);
  p0.c = (1.5*q, 1*q);
  p1.c = (2.5*q, 1*q);
  p2.c = (3.5*q, 1*q);
  p3.c = (4.5*q, 1*q);
  p4.c = (5.5*q, 1*q);
  p5.c = (6.5*q, 1*q);

  fill bpath.t1 withcolor .85white;
  fill bpath.t2 withcolor .85white;
  fill bpath.r1 withcolor .85white;
  fill bpath.r2 withcolor .85white;
  fill bpath.s1 withcolor .85white;

  drawboxed(u1, u2, i1, i2, i3, t1, t2, r1);
  drawboxed(r2, s1, p0, p1, p2, p3, p4, p5);
  
  draw u1.c--t1.c cutbefore bpath.u1 cutafter bpath.t1;
  draw u1.c--r1.c cutbefore bpath.u1 cutafter bpath.r1;
  draw t1.c--r2.c cutbefore bpath.t1 cutafter bpath.r2;
  draw r1.c--t2.c cutbefore bpath.r1 cutafter bpath.t2;
  draw r2.c--s1.c cutbefore bpath.r2 cutafter bpath.s1;
  draw t2.c--u2.c cutbefore bpath.t2 cutafter bpath.u2;
  draw s1.c--i1.c cutbefore bpath.s1 cutafter bpath.i1;
  draw u2.c--i2.c cutbefore bpath.u2 cutafter bpath.i2;
  draw u2.c--i3.c cutbefore bpath.u2 cutafter bpath.i3;
  draw i1.c--p0.c cutbefore bpath.i1 cutafter bpath.p0;
  draw i1.c--p1.c cutbefore bpath.i1 cutafter bpath.p1;
  draw i2.c--p2.c cutbefore bpath.i2 cutafter bpath.p2;
  draw i2.c--p3.c cutbefore bpath.i2 cutafter bpath.p3;
  draw i3.c--p4.c cutbefore bpath.i3 cutafter bpath.p4;
  draw i3.c--p5.c cutbefore bpath.i3 cutafter bpath.p5;
  
  drawarrow u1.c{left}..i1.c cutbefore bpath.u1 cutafter bpath.i1 dashed evenly;
  drawarrow u1.c{down}..{down}u2.c cutbefore bpath.u1 cutafter bpath.u2 dashed evenly;
endfig;

end
