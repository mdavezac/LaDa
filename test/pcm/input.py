""" Input file declaring the Point Charge + van der Wall model. """


clj  = Clj()
""" Point charge + vand der Wall model. """
clj.mesh = (5, 5, 5)
clj.ewald_cutoff = 25
clj.lj_cutoff = 10.5
clj.charges["Li"] = 1.0
clj.charges["Cs"] = 1.0
clj.charges["Br"] = -1.0
hs = { "Li":1.34, "Cs":2.25,"Br":1.14}
vdw = {"Li":2.2,"Cs":3,"Br":1.9}
for a in ["Li", "Cs", "Br" ]:
  for b in ["Li", "Cs", "Br" ]:
    type = models.bond_type( a, b )
    hs_ = float( hs[a] ) + float( hs[b]  )
    vdw_ = float(vdw[a]) + float(vdw[b] )
    clj.bonds[type] = models.LJBond( pow(hs_, 12.0), pow(vdw_, 6.0) )
