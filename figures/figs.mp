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

end
    