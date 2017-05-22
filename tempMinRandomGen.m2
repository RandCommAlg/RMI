randomMinGeneratingSet = method()
randomMinGeneratingSet(ZZ, ZZ, RR) := (n,D,p) -> (
    return randomMinGeneratingSet(n,D,D:p)
 )

randomMinGeneratingSet(ZZ, ZZ, Sequence) := (n,D,p) -> (
 x:=symbol x;
 R:=QQ[x_1..x_n];
 B:={};
 currentRing=R;
 apply(D, d-> (
  chosen:=select(flatten entries basis(d+1,currentRing), m->random(0.0,1.0)<=p_d);
  B=flatten append(B, chosen/(i->sub(i, R)));
  currentRing = currentRing/promote(ideal(chosen),currentRing)
  )
);
 return(B)
 )