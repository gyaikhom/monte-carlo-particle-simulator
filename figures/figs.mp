\input boxes

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
  label.rt(btex $index\_lookup\_table$ etex, (11.5, 4.5)*u);
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
  boxit.p0(btex 0 etex);
  boxit.p1(btex 1 etex);
  boxit.p2(btex 2 etex);
  boxit.p3(btex 3 etex);
  boxit.p4(btex 4 etex);
  boxit.p5(btex 5 etex);
  
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

outputtemplate := "geometry-solids.mps";
beginfig(1);
  numeric q;
  q := 0.75cm;
  circleit.u1(btex $\cup$ etex);
  circleit.u2(btex $\cup$ etex);
  circleit.u3(btex $\cup$ etex);
  circleit.i1(btex $\cap$ etex);
  circleit.i2(btex $\cap$ etex);
  circleit.d1(btex $\setminus$ etex);
  circleit.d2(btex $\setminus$ etex);

  boxit.p0(btex 0 etex);
  boxit.p1(btex 1 etex);
  boxit.p2(btex 2 etex);
  boxit.p3(btex 3 etex);
  boxit.p4(btex 4 etex);
  boxit.p5(btex 5 etex);
  boxit.p6(btex 6 etex);
  boxit.p7(btex 7 etex);
  boxit.p8(btex 8 etex);
  boxit.p9(btex 9 etex);

  p0.c = (1*q, q);
  p1.c = (2*q, q);
  p2.c = (3*q, q);
  p3.c = (4*q, q);
  u1.c = (1.5*q, 2*q);
  u2.c = (3.5*q, 2*q);
  i1.c = (2.5*q, 3*q);
  label.urt(btex $A$ etex, i1.c+(0.15*q,0.15*q));
  
  p4.c = (5*q, 2*q);
  p5.c = (6*q, 2*q);
  d1.c = (5.5*q, 3*q);
  label.urt(btex $B$ etex, d1.c+(0.2*q,0.2*q));

  p6.c = (7*q, q);
  p7.c = (8*q, q);
  p8.c = (9*q, q);
  p9.c = (10*q, q);
  u3.c = (7.5*q, 2*q);
  d2.c = (9.5*q, 2*q);
  i2.c = (8.5*q, 3*q);
  label.urt(btex $C$ etex, i2.c+(0.15*q,0.15*q));

  drawboxed(u1, u2, u3, i1, i2, d1, d2);
  drawboxed(p0, p1, p2, p3, p4, p5, p6, p7, p8, p9);

  draw i1.c--u1.c cutbefore bpath.i1 cutafter bpath.u1;
  draw i1.c--u2.c cutbefore bpath.i1 cutafter bpath.u2;
  draw u1.c--p0.c cutbefore bpath.u1 cutafter bpath.p0;
  draw u1.c--p1.c cutbefore bpath.u1 cutafter bpath.p1;
  draw u2.c--p2.c cutbefore bpath.u2 cutafter bpath.p2;
  draw u2.c--p3.c cutbefore bpath.u2 cutafter bpath.p3;

  draw d1.c--p4.c cutbefore bpath.d1 cutafter bpath.p4;
  draw d1.c--p5.c cutbefore bpath.d1 cutafter bpath.p5;

  draw i2.c--u3.c cutbefore bpath.i2 cutafter bpath.u3;
  draw i2.c--d2.c cutbefore bpath.i2 cutafter bpath.d2;
  draw u3.c--p6.c cutbefore bpath.u3 cutafter bpath.p6;
  draw u3.c--p7.c cutbefore bpath.u3 cutafter bpath.p7;
  draw d2.c--p8.c cutbefore bpath.d2 cutafter bpath.p8;
  draw d2.c--p9.c cutbefore bpath.d2 cutafter bpath.p9;

  draw (12,0)*q--(16,0)*q--(16,4)*q--(12,4)*q--cycle;
  draw (12,2)*q--(16,2)*q dashed evenly;
  draw (14,0)*q--(14,4)*q dashed evenly;

  path c[];
  c1 = fullcircle scaled 1q;
  c2 = fullcircle scaled 3q;
  c3 = fullcircle scaled 1q;
  draw c1 shifted (12.5*q,0.5*q);
  draw c2 shifted (14*q,2*q);
  draw c3 shifted (15.5*q,3.5*q);
  label(btex $A$ etex, (12.5*q,0.5*q));
  label(btex $B$ etex, (13.5*q,2.5*q));
  label(btex $C$ etex, (15.5*q,3.5*q));

  label.lrt("0", (12*q,4*q));
  label.lrt("1", (14*q,4*q));
  label.lrt("2", (12*q,2*q));
  label.lrt("3", (14*q,2*q));
endfig;

outputtemplate := "geometry-primitives-table.mps";
beginfig(1);
  numeric q;
  q := 0.75cm;
  label(btex Primitives Table etex, (2,5.75)*q);
  draw (0,0)*q--(4,0)*q--(4,6)*q--(0,6)*q--cycle;
  draw (0,5.5)*q--(4,5.5)*q;
  draw (1.25,5.5)*q--(1.25,0)*q;
  draw (2.75,5.5)*q--(2.75,0)*q;
  label(btex $a$ etex, (0.5[0,1.25],5.25)*q);
  label(btex $i$ etex, (0.5[1.25,2.75],5.25)*q);
  label(btex $p$ etex, (0.5[2.75,4],5.25)*q);
  numeric j, p[];
  j := 0.5;
  for i=0 upto 9:
    draw (0,j)*q--(4,j)*q;
    p[i] := j-0.25;
    label.lft(decimal(9-i), (0,p[i])*q);
    j := j + 0.5;
  endfor
  draw (21.25,5.5)*q--(5,5.5)*q--(5,5)*q--(21.25,5)*q--cycle;
  numeric v[];
  j := 8.5;
  for i=0 upto 16:
    draw (j,5.5)*q--(j,5)*q;
    v[i] := j + 0.375;
    label(decimal(i), (v[i],5.75)*q);
    j := j + 0.75;
  endfor
  label(btex Postfix Buffer etex, (0.5[5,8.5],5.25)*q);
  label("0", (v[0],5.25)*q);
  label("1", (v[1],5.25)*q);
  label("-3", (v[2],5.25)*q);
  label("2", (v[3],5.25)*q);
  label("3", (v[4],5.25)*q);
  label("-3", (v[5],5.25)*q);
  label("-2", (v[6],5.25)*q);
  label("4", (v[7],5.25)*q);
  label("5", (v[8],5.25)*q);
  label("-1", (v[9],5.25)*q);
  label("6", (v[10],5.25)*q);
  label("7", (v[11],5.25)*q);
  label("-3", (v[12],5.25)*q);
  label("8", (v[13],5.25)*q);
  label("9", (v[14],5.25)*q);
  label("-1", (v[15],5.25)*q);
  label("-2", (v[16],5.25)*q);
  label(btex $\cup$ etex, (v[2],4.75)*q);
  label(btex $\cup$ etex, (v[5],4.75)*q);
  label(btex $\cap$ etex, (v[6],4.75)*q);
  label(btex $\setminus$ etex, (v[9],4.7)*q);
  label(btex $\cup$ etex, (v[12],4.75)*q);
  label(btex $\setminus$ etex, (v[15],4.7)*q);
  label(btex $\cap$ etex, (v[16],4.75)*q);
  drawarrow (v[0],5)*q{down}..tension 2.5..{left}(4,p[9])*q;
  drawarrow (v[1],5)*q{down}..tension 2.5..{left}(4,p[8])*q;
  drawarrow (v[3],5)*q{down}..tension 2.5..{left}(4,p[7])*q;
  drawarrow (v[4],5)*q{down}..tension 2.5..{left}(4,p[6])*q;
  drawarrow (v[7],5)*q{down}..tension 2.5..{left}(4,p[5])*q;
  drawarrow (v[8],5)*q{down}..tension 2.5..{left}(4,p[4])*q;
  drawarrow (v[10],5)*q{down}..tension 2.5..{left}(4,p[3])*q;
  drawarrow (v[11],5)*q{down}..tension 2.5..{left}(4,p[2])*q;
  drawarrow (v[13],5)*q{down}..tension 2.5..{left}(4,p[1])*q;
  drawarrow (v[14],5)*q{down}..tension 2.5..{left}(4,p[0])*q;
endfig;

outputtemplate := "geometry-solids-table.mps";
beginfig(1);
  numeric q;
  q := 0.75cm;
  label(btex Solids Table etex, (3.75,5.75)*q);
  draw (1.5,3.5)*q--(6,3.5)*q--(6,6)*q--(1.5,6)*q--cycle;
  draw (1.5,5.5)*q--(6,5.5)*q;
  draw (2.75,5.5)*q--(2.75,3.5)*q;
  draw (4.375,5.5)*q--(4.375,3.5)*q;
  label(btex $bb$ etex, (0.5[1.5,2.75],5.25)*q);
  label(btex $start$ etex, (0.5[2.75,4.375],5.25)*q);
  label(btex $count$ etex, (0.5[4.375,6],5.25)*q);
  numeric j, p[];
  j := 4;
  for i=0 upto 2:
    draw (1.5,j)*q--(6,j)*q;
    p[i] := j-0.25;
    label.lft(decimal(2-i), (1.5,p[i])*q);
    j := j + 0.5;
  endfor
  draw (23.25,5.5)*q--(7,5.5)*q--(7,5)*q--(23.25,5)*q--cycle;
  numeric v[];
  j := 10.5;
  for i=0 upto 16:
    draw (j,5.5)*q--(j,5)*q;
    v[i] := j + 0.375;
    label(decimal(i), (v[i],5.75)*q);
    j := j + 0.75;
  endfor
  label(btex Postfix Buffer etex, (0.5[7,10.5],5.25)*q);
  label("0", (v[0],5.25)*q);
  label("1", (v[1],5.25)*q);
  label("-3", (v[2],5.25)*q);
  label("2", (v[3],5.25)*q);
  label("3", (v[4],5.25)*q);
  label("-3", (v[5],5.25)*q);
  label("-2", (v[6],5.25)*q);
  label("4", (v[7],5.25)*q);
  label("5", (v[8],5.25)*q);
  label("-1", (v[9],5.25)*q);
  label("6", (v[10],5.25)*q);
  label("7", (v[11],5.25)*q);
  label("-3", (v[12],5.25)*q);
  label("8", (v[13],5.25)*q);
  label("9", (v[14],5.25)*q);
  label("-1", (v[15],5.25)*q);
  label("-2", (v[16],5.25)*q);
  drawarrow (6,p[2])*q{right}..tension 3..{up}(v[0],5)*q;
  drawarrow (6,p[1])*q{right}..tension 3..{up}(v[7],5)*q;
  drawarrow (6,p[0])*q{right}..tension 3..{up}(v[10],5)*q;
  label("0", (0.5[2.75,4.375],p[2])*q);
  label("7", (0.5[4.375,6],p[2])*q);
  label("7", (0.5[2.75,4.375],p[1])*q);
  label("3", (0.5[4.375,6],p[1])*q);
  label("10", (0.5[2.75,4.375],p[0])*q);
  label("7", (0.5[4.375,6],p[0])*q);
  label(btex $A$ etex, (0.5[v[0],v[7]],4.7)*q);
  label(btex $B$ etex, (0.5[v[7],v[10]],4.7)*q);
  label(btex $C$ etex, (0.5[v[10],v[16]],4.7)*q);
endfig;

outputtemplate := "geometry-subcuboids-table.mps";
beginfig(1);
  numeric q;
  q := 0.75cm;
  label(btex Subcuboids Table etex, (3,5.75)*q);
  draw (0,3)*q--(6,3)*q--(6,6)*q--(0,6)*q--cycle;
  draw (0,5.5)*q--(6,5.5)*q;
  draw (2,5.5)*q--(2,3)*q;
  draw (4,5.5)*q--(4,3)*q;
  label(btex $bb$ etex, (1,5.25)*q);
  label(btex $start$ etex, (3,5.25)*q);
  label(btex $count$ etex, (5,5.25)*q);
  numeric j, p[];
  j := 3.5;
  for i=0 upto 3:
    draw (0,j)*q--(6,j)*q;
    p[i] := j-0.25;
    label.lft(decimal(3-i), (0,p[i])*q);
    j := j + 0.5;
  endfor
  draw (15,5.5)*q--(7,5.5)*q--(7,5)*q--(15,5)*q--cycle;
  numeric v[];
  j := 10.5;
  for i=0 upto 5:
    draw (j,5.5)*q--(j,5)*q;
    v[i] := j + 0.375;
    label(decimal(i), (v[i],5.75)*q);
    j := j + 0.75;
  endfor
  label(btex Solid Indices etex, (0.5[7,10.5],5.25)*q);
  label("1", (v[0],5.25)*q);
  label("1", (v[1],5.25)*q);
  label("2", (v[2],5.25)*q);
  label("0", (v[3],5.25)*q);
  label("1", (v[4],5.25)*q);
  label("1", (v[5],5.25)*q);
  drawarrow (6,p[3])*q{right}..tension 3..{up}(v[0],5)*q;
  drawarrow (6,p[2])*q{right}..tension 3..{up}(v[1],5)*q;
  drawarrow (6,p[1])*q{right}..tension 3..{up}(v[3],5)*q;
  drawarrow (6,p[0])*q{right}..tension 3..{up}(v[5],5)*q;
  label("0", (3,p[3])*q);
  label("1", (5,p[3])*q);
  label("1", (3,p[2])*q);
  label("2", (5,p[2])*q);
  label("3", (3,p[1])*q);
  label("2", (5,p[1])*q);
  label("5", (3,p[0])*q);
  label("1", (5,p[0])*q);
endfig;

outputtemplate := "node-type-bit-field.mps";
beginfig(1);
  numeric q;
  q := 0.45cm;
  draw (0,3)*q--(32,3)*q--(32,4)*q--(0,4)*q--cycle;
  draw (28,3)*q--(28,4)*q;
  draw (29,3)*q--(29,4)*q;
  label(btex Line number in geometry file which generated the node etex, (14,3.5)*q);
  label(btex Usage bit etex, (26.5,4.75)*q);
  label(btex Node type etex, (30.5,4.75)*q);
  drawarrow (26.5,4.25)*q--(28.5,3.5)*q;
  drawarrow (30.5,4.25)*q--(30.5,3.5)*q;
  label("31", (0.5,2.5)*q);
  label("3", (28.5,2.5)*q);
  label("0", (31.5,2.5)*q);
endfig;

end
