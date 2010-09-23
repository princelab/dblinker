

def putsv(*args)
	puts(*args) if $VERBOSE
end

# Narray based distance calculator
#~ module Bio
	#~ class PDB
		#~ module Utils
			def distance_checker(coord, other, cut_off)
			sq_diffs=[]
				(0...(coord.size)).each do |i|
				pos = coord[i]
				oth = other[i]
				sq_diffs << (oth - pos)**2
				end
			NMath.sqrt(sq_diffs.inject {|sum, vec| sum + vec }) #< cut_off
			end
		#~ end
	#~ end
#~ end
=begin
def dist(coord, xyz_arr, cut_off)
	# coord is an array of coordinate arrays
	# xyz_arr is an array with all possible xyz coordinates in it
	# cut_off is the defined cutoff for early exit checking
	# output will an array of trues corresponding to how many of coord were within the distance cutoff
	coord.each do |pt|
		sq_diffs=[pt.size]
		(0...(pt.size).each do |i|
			pos = coord[i]
			xyz_arr[i].each do |j|
			sq_diffs[i] << (oth-j)**2
end
end
end
end

=end			

# Quiet the output of the PDB read
def quiet(&block)
	old_level = $VERBOSE
	reply = block.call
	$VERBOSE = old_level
	reply
end
# Tunneler fxn
def tunneler (vec1, vec2, num_of_steps)
	result = []
	increment = [(vec2[0] - vec1[0])/(num_of_steps+1), (vec2[1]-vec1[1])/(num_of_steps+1), (vec2[2]-vec1[2])/(num_of_steps+1)]
	curr_point = [vec1[0]+increment[0], vec1[1]+increment[1],vec1[2]+increment[2]]
	i=0
	while i < num_of_steps
		curr_vec = []
		j=0
		while j < 3 
			curr_vec[j] = (curr_point[j])		
			curr_point[j] += increment[j]
			j += 1
		end
		result[i] = curr_vec
		i +=1
	end
	return result
end


#class Array
    # A binary search library adapted from: http://0xcc.net/ruby-bsearch/
    # ---
    #
    # Ruby/Bsearch - a binary search library for Ruby.
    #
    # Copyright (C) 2001 Satoru Takabayashi <satoru@namazu.org>
    #     All rights reserved.
    #     This is free software with ABSOLUTELY NO WARRANTY.
    #
    # You can redistribute it and/or modify it under the terms of 
    # the Ruby's licence.
    # 
    # The binary search algorithm is extracted from Jon Bentley's
    # Programming Pearls 2nd ed. p.93
    #
 
      #
      # Return the lower boundary. (inside)
      #
      def search_lower_boundary(array, range=nil, &block)
        range = 0 ... array.length if range == nil
        lower  = range.first() -1
        upper = if range.exclude_end? then range.last else range.last + 1 end
        while lower + 1 != upper
          mid = ((lower + upper) / 2).to_i # for working with mathn.rb (Rational)
          if yield(array[mid]) < 0
            lower = mid
          else 
            upper = mid
          end
        end
        return upper
      end
    
      #
      # Return the upper boundary. (outside)
      #
      def search_upper_boundary(array, range = nil, &block)
        range = 0 ... array.length if range == nil        
        lower  = range.first() -1
        upper = if range.exclude_end? then range.last else range.last + 1 end
        while lower + 1 != upper
          mid = ((lower + upper) / 2).to_i # for working with mathn.rb (Rational)
          if yield(array[mid]) <= 0
            lower = mid
          else 
            upper = mid
          end
        end
        return lower + 1 # outside of the matching range.
      end
#end

	#
	# Ryan's fxn for finding the window
	#
		def bin_window(array, point, int, cutoff)
			range = nil; 
			lower = search_lower_boundary(array) {|x| x[int] <=> (point - cutoff) } 
			upper = array.search_upper_boundary(array) {|x| x[int] <=> (point + cutoff) }
			return lower ... upper
		end
	   
	#
	# Ryan's fxn for finding the box 
	#
		def bin_box(array, min, max, int, cutoff)
			lower = search_lower_boundary(array) {|x| x[int] <=> (min - cutoff) } 
			upper = array.search_upper_boundary(array) {|x| x[int] <=> (max + cutoff) }
			return lower ... upper
		end  

# Linker class
class Linker < Array	
	# takes (atom1, atom2, accession, inter?, inside, steps, dist)
	def initialize(*args)
		super(9)
		self.replace(args)
		inside = nil
	end
	def inter?	
		if self[3].nil?
			self[3] = (self[0].chainID != self[1].chainID)
		else
			self[3]
		end
	end
	def file
		if self[2].nil? 		# Str.accession
			self[2] = 'undefined'
		else
			self[2]
		end
	end
	def atom1
		self[0]
	end
	def atom2 
		self[1]
	end
	def inside
		if self[4].nil?
			self[4] = 'undefined'
		else
			self[4]
		end
	end
	def steps
		if self[5].nil?
			self[5] = 'undefined'
		else
			self[5]
		end
	end
	def x
		[atom1.x, atom2.x]
	end
	def y
		[atom1.y, atom2.y]
	end
	def z
		[atom1.z, atom2.z]
	end
	def vec1
		atom1.xyz
	end
	def vec2
		atom2.xyz
	end
	def out
		#~ puts 'Atom chains are: '+  atom1.chainID  + ' and '+ atom2.chainID
		puts 'Residues are: ' + atom1.resName + ' ' + atom1.resSeq.to_s + ' and ' + atom2.resName + ' ' + atom2.resSeq.to_s
		puts 'Atom coordinates are: ' + self.vec1.to_s + ' and ' + self.vec2.to_s
		puts 'Distance is: ' + self.dist.to_s
		puts '% insideness is: ' + inside.to_s; 
	end
	def dist
		if self[6].nil?
			self[6]=vec1.distance(vec2)
		else
			self[6]
		end
	end
	def solv_access_1	# this represents the minimum dist to the solvent
		self[7]
	end
	def solv_access_2
		self[8]
	end
	def to_archive(*args)
		arch = { "PDB" => file, "inter" => inter?, 
			"insideness: " => "#{inside*100} %"
		}
		arch.to_yaml(*args)
	end
	def to_yaml(*args)
		hash = { "ID_1" => atom1.serial, "ID_2" => atom2.serial,
			"chain1" => atom1.chainID, "chain2" => atom2.chainID,
			"coord1" => vec1.to_a, "coord2" => vec2.to_a,
			"dist" => dist,
			"inter" => inter? ,
			"inside" => inside,
			"steps" => steps,
			"dist_to_solvent_1" => solv_access_1,
			"dist_to_solvent_2" => solv_access_2
		}
		hash.to_yaml(*args)
	end
end
	
