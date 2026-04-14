# ligands

mol representation CPK
mol selection {resname ABU Y4B}
mol color Element
mol material AOChalky
mol addrep top

# protein chains

mol representation NewCartoon
mol selection {protein and chain A}
mol color ColorID 31
mol material AOChalky
mol addrep top

mol representation NewCartoon
mol selection {protein and chain B}
mol color ColorID 23
mol material AOChalky
mol addrep top

mol representation NewCartoon
mol selection {protein and chain C}
mol color ColorID 9
mol material AOChalky
mol addrep top

mol representation NewCartoon
mol selection {protein and chain D}
mol color ColorID 15
mol material AOChalky
mol addrep top

mol representation NewCartoon
mol selection {protein and chain E}
mol color ColorID 2
mol material AOChalky
mol addrep top

# protein residues

mol representation CPK
mol selection {chain B D and resid 46 67 and not element H}
mol color Element
mol material AOChalky
mol addrep top

mol representation CPK
mol selection {chain A C and resid 155 200 and not element H}
mol color Element
mol material AOChalky
mol addrep top

# selected glycans

mol representation Licorice
mol selection {glycan or (chain A C and resid 149 80)}
mol color ColorID 2
mol material AOChalky
mol addrep top

