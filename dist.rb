#!/usr/bin/env ruby 

# Develop dist code for insertion into xlink_dist_calc.rb (Final, hopefully, code)




	# John's Code

	  pdb_with_hydrogens =
	    if opt[:add_hydrogen]
	      pdb_plus_h_added = base + '_Hadded.pdb'
	      putsv "writing to: #{pdb_plus_h_added}"
	      Pymol::HydrogenBonds.pdb_with_hydrogens(file, pdb_plus_h_added)
	    else
	      file
	    end

	  base_h_added = pdb_with_hydrogens.chomp(File.extname(pdb_with_hydrogens))

	  hbond_arrays = Pymol::HydrogenBonds.from_pdb(pdb_with_hydrogens, opt)

	  # http://pymolwiki.org/index.php/Surface#Exporting_Surface.2FMesh_Coordinates_to_File
	  surface_coords = Pymol::Surface.from_pdb(pdb_with_hydrogens)

	  sc_sz = surface_coords.size
	  (xs, ys, zs) = [nil,nil,nil].map { NArray.float(sc_sz) }
	  surface_coords.each_with_index do |xyz, i|
	    xs[i] = xyz[0]
	    ys[i] = xyz[1]
	    zs[i] = xyz[2]
	  end

	  # get the distance from hydrogen to surface
	  # 0 => donor
	  # 1 => hydrogen
	  # 2 => acceptor
	  which_atom = 1

	  # just output distance to the surface of the amino acid
	  if opt[:aa_to_surf]
	    amino_acids = []
	    hbond_arrays.each do |a,b,c|
	      [a,c].each {|atom| amino_acids << atom.residue }
	    end
	    amino_acids.uniq.each do |res|

	      #### CENTER OF GRAVITY
	      #coord = res.centreOfGravity
	      #### GEOMETRIC CENTER
	      #coord = res.geometricCentre
	      #min_dist = Bio::PDB::Utils.distance_to_many(coord, [xs, ys, zs] ).min

	      #### MINIMUM DISTANCE TO ANY ATOM
	      min_dist = res.atoms.map {|atom| Bio::PDB::Utils.distance_to_many(atom.xyz, [xs, ys, zs] ).min }.min
	      # output to mimic older output
	      puts "#{res.resName}#{res.id} minimum_distance_to_surface #{min_dist}"
	    end

	    next
	  end

	  putsv "calculating distances to surface ..."
	  # also we are gathering all the data we need.
	  characterized = hbond_arrays.map do |data|
	    coords = Array.new(3)
	    na_coords = Array.new(3)
	    data[0,3].each_with_index do |atom,i|
	      coords[i] = atom.xyz
	      na_coords[i] = NArray.to_na(coords[i].to_a)
	    end

	    dists_to_surface = Bio::PDB::Utils.distance_to_many(coords[which_atom], [xs, ys, zs] )

	    data[3] = Bio::PDB::Utils.rad2deg(data[3]) unless opt[:radians]
	    ids = data[0,3].map {|atom| atom.serial }
	    (don, acc) = [data[0], data[2]].map {|at| [at.resName, at.residue.id, at.name] }
	    ids.push(data[1].name)
	    id_part = ids.push(*don).push(*acc)
	    id_part.push(*(data[3,3]))
	    id_part.push(dists_to_surface.min)
	  end

	  final_output = base_h_added + output_postfix
	  File.open(final_output, 'w') do |out|
	    out.puts categories.join(opt[:delim])
	    characterized.each do |array|
	      out.puts array.join(opt[:delim])
	    end
	  end
