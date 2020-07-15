
def number_permutation_groups_up_to_degree(N):
  # Counting will be performed for degrees (domain size) from 1 to N+1
  # Returns the sequence counting permutation groups (up to isomorphism) as a list
  res = []
  for i in range(N+1):  # initialisation, with the trivial groups on all degrees
    res.append(1)
  for total_deg in range(2, N+2):
    for pt in list(Partitions(total_deg)):  # set partition into orbits (up to conj.)
      l = len(pt)
      pool_of_gps = [] # eventually the set of possibilities of restrictions to the orbits
      i = 0
      triv_orb = 0     # number of orbits that are singletons (``ignored'')
      ambient_gp = PermutationGroup([]) # will be used to check the conjugacy later
      while i < l:
        #### A. Creation of a set of vectors that represent all possible direct products on these orbits
        # The coordinates will be indices of groups in libgap.TransitiveGroups
        prec = pt[i]                  # size of the preceding part/orbit
        if prec == 1:                 # the group is trivial on remaining domain
          triv_orb = l-i              # number of trivial orbits, at the end of pt
          break                       # only one choice, no need to go further
        nb_rep = 0                    # nb of orbits of same size (repetitions)
        while i < l and pt[i] == prec:
          # identifies same size orbits to avoid isomorphic choices of groups (eg. G1*G2 and G2*G1)
          nb_rep += 1
          i += 1
        factor = libgap.WreathProduct(SymmetricGroup(prec), SymmetricGroup(nb_rep))
        ambient_gp = libgap.DirectProduct(ambient_gp, factor)
        Ntrans = libgap.NrTransitiveGroups(prec)
        S = SymmetricGroup(nb_rep)
        pool_of_gps.append(IntegerVectorsModPermutationGroup(S, max_part=Ntrans-1))
      pool_of_gps = cartesian_product(pool_of_gps)  # set of (tuples of) vectors
      #### B. Computing of all possible direct products with the given orbits (pt)
      if triv_orb == l:           # trivial partition so trivial group
        continue                  # continue to next partition
      for elt in pool_of_gps:
        # each elt represents a direct product of transitive permutation groups
        gps = []   # the actual factors of the product (isom. to the orbit restrictions)
        ind_orb = 0
        # creation of the direct product (gap element) from which to consider the subdirect products
        for vector in elt:
          # each vector represents the choices of groups on a given set of same-size orbits
          for n in vector:
            gps.append(libgap.TransitiveGroup(pt[ind_orb], n+1))
            ind_orb += 1
        #### C. Computing of all subdirect products of the groups in gps (ignores trivial ones), up to conjugacy
        subs = [gps[0]];                  # subdirect products of j groups from gps
        for j in range(1, l-triv_orb):    # l-triv_orb = len(gps)
          new_subs = []            # will include one more group in the subdirect products
          for sub in subs:
            candidates = list(libgap.SubdirectProducts(sub, gps[j]))
            for new in candidates:        # checking conjugacy
              conjugate = false
              for other in new_subs:
                if libgap.IsConjugate(ambient_gp, new, other):
                  conjugate = true
                  break
              if not conjugate:
                new_subs.append(new)
          subs = new_subs
        res[total_deg-1] += len(subs)
  return res
    
