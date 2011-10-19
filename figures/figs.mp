vardef arrowhead expr p =
  save A,u,a,b; pair A,u; path a,b;
  A := point length(p) of p;
  u := unitvector(direction length(p) of p);
  a := A{-u} .. (A - ahlength*u rotated 20);
  b := A{-u} .. (A - ahlength*u rotated -20);
  ( a & reverse(a) & b & reverse(b) ) --cycle
enddef;

outputtemplate := "neighbour-direct.mps";
beginfig(1);
  numeric u;
  u = 0.5cm;
  label.bot(btex Direct Mapping etex, (0.5[13,0], 0)*u);
  label.lft("1", (0, 1.5)*u);
  label.lft("0", (0, 2.5)*u);
  label.rt(btex $s$-field etex, (13.5, 4.5)*u);
  label.rt(btex $neighbour\_table$ etex, (13.5, 2.5)*u);
  fill (2,4)*u--(3,4)*u--(3,5)*u--(2,5)*u--cycle withcolor .85white;
  fill (6,4)*u--(7,4)*u--(7,5)*u--(6,5)*u--cycle withcolor .85white;
  draw (0,0)*u--(0,3)*u--(13,3)*u--(13,0)*u;
  for r=1 upto 3: draw (0,r)*u--(13,r)*u; endfor
  draw (0,4)*u--(13,4)*u--(13,5)*u--(0,5)*u--cycle;
  for c=1 upto 8:
    draw (c,4)*u--(c,5)*u;
    label(decimal(c),(c-0.5,4.5)*u);
  endfor
  label(btex $\cdots$ etex,(9.5,4.5)*u);
  for c=11 upto 12:
    draw (c,4)*u--(c,5)*u;
    label(decimal(c + 30),(c+0.5,4.5)*u);
  endfor
  for c=0.5 upto 1.5:
    drawarrow (c,4)*u--(c,3)*u;
  endfor
  for c=3.5 upto 6:
    drawarrow (c,4)*u--(c,3)*u;
  endfor
  drawarrow (7.5,4)*u--(7.5,3)*u;
  for c=11.5 upto 12.5:
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
  label.rt(btex $s$-field etex, (13.5, 6.5)*u);
  label.rt(btex $neighbour\_idx\_table$ etex, (11.5, 4.5)*u);
  label.rt(btex $neighbour\_table$ etex, (11.5, 2.5)*u);
  draw (0,4)*u--(11,4)*u--(11,5)*u--(0,5)*u--cycle;
  for c=1 upto 6:
    draw (c,5)*u--(c,4)*u;
    drawarrow (c-0.5,4)*u--(c-0.5,3)*u;
    label(decimal(c-1),(c-0.5,4.5)*u);
  endfor
  label(btex $\cdots$ etex,(7.5,4.5)*u);
  for c=9 upto 10:
    draw (c,5)*u--(c,4)*u;
    drawarrow (c+0.5,4)*u--(c+0.5,3)*u;
    label(decimal(c+15),(c+0.5,4.5)*u);
  endfor
  draw (0,0)*u--(0,3)*u--(11,3)*u--(11,0)*u;
  for r=1 upto 3: draw (0,r)*u--(11,r)*u; endfor

  fill (2,6)*u--(3,6)*u--(3,7)*u--(2,7)*u--cycle withcolor .85white;
  fill (6,6)*u--(7,6)*u--(7,7)*u--(6,7)*u--cycle withcolor .85white;
  draw (0,6)*u--(13,6)*u--(13,7)*u--(0,7)*u--cycle;
  for c=1 upto 8:
    draw (c,6)*u--(c,7)*u;
    label(decimal(c),(c-0.5,6.5)*u);
  endfor
  label(btex $\cdots$ etex,(9.5,6.5)*u);
  for c=11 upto 12:
    draw (c,6)*u--(c,7)*u;
    label(decimal(c + 30),(c+0.5,6.5)*u);
  endfor
  for c=0.5 upto 1.5:
    drawarrow (c,6)*u--(c,5)*u;
  endfor
  for c=3.5 upto 6:
    drawarrow (c,6)*u--(c-1,5)*u;
  endfor
  drawarrow (7.5,6)*u--(5.5,5)*u;
  for c=11.5 upto 12.5:
    drawarrow (c,6)*u--(c-2,5)*u;
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

end
    