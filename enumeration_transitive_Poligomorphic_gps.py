
def enumeration_transitive_Poligomorphic_gps(N) :
  res = []
  for i in range(N+1) :  # initialisation; counting will be performed for growths from 0 to N
    res.append(0)
  # trivial case: degree 1, for libgap.TransitiveGroup needs the degree to be > 1
  res[0] += 5  # hands the 5 highly homogeneous groups
  for deg in range(2, N+2) :  # a tight lower bound on the growth is deg - 1
    for i in range(libgap.NrTransitiveGroups(deg)) :
      gp = libgap.TransitiveGroup(deg, i+1)
      # treatment of the two trivial block systems : first, blocks of size 1
      res[deg-1] += 5  # 5 choices of highly homogeneous group in that case (and growth = deg-1)
      # second, one single block (will hand a single superblock)
      for H in libgap.NormalSubgroups(gp) :
        H_sage = PermutationGroup_generic(gap_group = H)
        if H_sage == PermutationGroup_generic([]) : # for else sage assumes the domain to be {1}
          H_sage = PermutationGroup_generic([], domain = range(1, deg+1))
        if len(H_sage.domain()) <> deg:
          print("Warning: domain likely narrowed during the conversion from gap to sage.")
        series_H = H_sage.profile_series()
        nb_orbits_H = sum(series_H[n] for n in range(1, deg+1))  # total nb of non trivial orbits of H
        growth = nb_orbits_H - 1  # growth = nb of gens - 1 = non trivial orbits - 1
        if growth <= N :        
          res[growth] += 1
      # dom = H_sage.domain()   # {1, 2, ... , deg} 
      # non trivial block systems
      for repres in libgap.AllBlocks(gp) :  # a representative of blocks for each system
        stab = libgap.Stabilizer(gp, repres, libgap.OnSets)
        G = libgap.Image(libgap.ActionHomomorphism(stab, repres))  # restriction to the block
        size = libgap.Size(repres)
        for H in libgap.NormalSubgroups(G) :
          H_sage = PermutationGroup_generic(gap_group = H)
          if H_sage == PermutationGroup_generic([]) : # for else sage assumes the domain to be {1}
            H_sage = PermutationGroup_generic([], domain = range(1, size+1))
          if len(H_sage.domain()) <> size:
            print("Warning: domain likely narrowed during the conversion from gap to sage.")
          series_H = H_sage.profile_series()
          nb_orbits_H = sum(series_H[n] for n in range(1, size+1))  # non trivial orbits of H
          growth = nb_orbits_H*(deg / size) - 1  # non trivial orbits * nb of blocks - 1 = nb of gens - 1
          if growth <= N :
            res[growth] += 1
  return res

# 5, 6, 14, 33, 32, 114, 47, 323, 260, 338, 50, 2108, 58, 430, 940, 12470, 60, 7361, 64, 12136
