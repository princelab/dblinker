#!/usr/bin/ruby 

# SUPERCOMPUTER STUFF
# #!/fslgroups/fslg_princelab/programs/bin/ruby
# Folder variable ID #
# Num=Dir.pwd[/[0-9]$/]
# exit if Num.size != 1

# Gem requirements
require 'rubygems'
require 'bio'
require 'yaml'
require 'narray'
require 'pymol'
require 'lib/pymol/surface'
require 'lib/pymol/connections'
require 'lib/pymol/hydrogen_bonds'
require 'lib/hydrogen_bondifier/utils'
require 'bio/db/pdb'

# Local (User) library
require 'ryan_lib'
require 'binary_search'

# Constants
Cut_off = 4.0
Step_dist = 3.0
Outer_cutoff = 30.0
Num=1

# Directory variable
Out_dir="o_#{Num}"

# Usage guidelines
if ARGV.size == 0
	puts "usage: #{File.basename(__FILE__)} < file>.pdb ... "
	puts "output: <file>.allowed_xlinks.yml ... "
	exit
end
#~  Essentially the xlink_finder program...
ARGV.each do |file|
# Create the *.working file
 a = File.exist?("#{file.sub(/\.pdb$/i, ".working")}")
 puts "#{file} has already been opened for analysis" if a
 system "touch #{file.sub(/\.pdb$/i, ".working")}"
 next if a 
# Check if file is already done...
# b = File.exist?("../#{Out_dir}/#{file.sub(/\.pdb$/i, ".xlink.yml")}")
 b = File.exist?("#{file.sub(/\.pdb$/i, ".xlink.yml")}")
 puts "#{file} has already been x-linked" if b
 next if b

# Surface Generation
surface_coords = Pymol::Surface.from_pdb(file)
	  sc_sz = surface_coords.size
	  (xs, ys, zs) = [nil,nil,nil].map { NArray.float(sc_sz) }
	  surface_coords.each_with_index do |xyz, i|
	    xs[i] = xyz[0]
	    ys[i] = xyz[1]
	    zs[i] = xyz[2]
	  end

#PDB Parser
	str = quiet do 
		Bio::PDB.new(IO.read(file))
	end
	amine_nitrogens = str.find_atom {|a| a.name.upcase == 'NZ'}
	linkers = []
	amine_nitrogens.each do |atom1|
	amine_nitrogens.each do |atom2|
		next if atom1 == atom2
		link =  Linker.new(atom1,atom2, str.accession)
		b = link.dist 
		next if b > Outer_cutoff || b < Cut_off
		link[5] = (b/Step_dist).ceil # 5 is the index for steps inside the Linker	
		link[7] = Bio::PDB::Utils.distance_to_many(atom1.xyz, [xs,ys,zs] ).min
		link[8] = Bio::PDB::Utils.distance_to_many(atom2.xyz, [xs,ys,zs] ).min
		linkers << link
	end
	end
#~ Dist_to_many fxn prep (building lists of all atoms)
	coordx = []; allx=[];ally=[];allz=[];
#~ atom_list = Time.now
	str.each_atom do |atom|
		allx << atom.x
		ally << atom.y
		allz << atom.z
	end
#~ puts "atom list generation took #{Time.now-atom_list}"
	all = [allx,ally,allz]
#~ Prep for indexed sorting for union
	x_set=[];# y_set=[]; z_set=[];
	i=0;	while i < allx.size
		x_set << [allx[i],i]
		i +=1
	end
# Sorting
	x_set.sort!; #y_set.sort!; z_set.sort!	;
# Distance checker
	i=0
	# Iterates through each linker
	while i < linkers.size
		steps = linkers[i].steps.to_f
	# For each step through vector space along the link, determines the distance to all other atoms
	dist = tunneler(linkers[i].vec1.to_a, linkers[i].vec2.to_a, steps).map do |j| 
	# X values
		out = x_set[bin_window(x_set, j[0], 0, Cut_off)].sort_by{|x| x[1]}.map{|a| a[1]}
	# Early Exit	
		if out.empty? # exits early if there aren't any hits
			false
		else
		# Generation of pre-sorted coords for dist_check
			coords = all.map do |arr|  
				out.map do |index| 
					arr[index] 
				end
			end 
		# Make everything NArray
			coords.map! do |arr|
				NArray.to_na(arr)
			end 
		# Run distance calculation on them
			(distance_checker(j, coords, Cut_off) < Cut_off).any?
		end
	end
	# Cleans up the dist array to measure the number of true entries
		dist.delete(false); 
	# Writes out the information to the parent linker
		linkers[i][4] = (dist.size/steps)
	# Increment
		i += 1
	end								
# Linkers output to yaml	
	ymlfile1 = file.sub(/\.pdb$/i, ".xlink.yml")
	puts "creating: #{ymlfile1}" if $VERBOSE
	File.open("#{ymlfile1}", "w") {|out| out.print linkers.sort.to_yaml }
# ####File.open("../#{Out_dir}/#{ymlfile1}", "w") {|out| out.print linkers.sort.to_yaml }
# Inside added to database file
	#~ out2file = "outs.yml"
	#~ short_out = linkers.map do |a| a.to_archive end 
	#~ File.open(out2file, "w") {|output| output.print short_out.to_yaml }
system "rm #{file.sub(/\.pdb$/i, ".working")}"
end
