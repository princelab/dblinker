#!/usr/bin/env ruby

require 'yaml'
Cut_off=0.8
Cylinder_diameter=0.1
# Usage guidelines
if ARGV.size == 0
	puts "Usage: #{File.basename(__FILE__)} <file>.xlink.yml ... "
	puts 'Output: <file>.xlink.yml.png & <file>.xlink.yml.pse'
	exit
end
Link=Struct.new(:atom1, :atom2)
ARGV.each do |file|

	input = YAML::load(IO.read(file))
	outs = input.map do |link|
		atom1= link.fetch("ID_1")
		atom2= link.fetch("ID_2")
		if link.fetch("inside") < Cut_off
			Link.new(atom1, atom2)
		else
			nil
		end
	end
outs.compact!
	filename = File.basename(file, ".xlink.yml")
	%x[rm tmp.pml] if File.exists?('tmp.pml')
	tf = File.open('tmp.pml','w+')
#	tf.puts 'log_open ("tmp_log.txt")'
	tf.puts 'run draw_links.py'
	tf.puts "load #{filename}.pdb"
	tf.puts 'hide everything'
	tf.puts 'show cartoon'
outs.each do |a|
puts a.atom2
puts a.atom1
	tf.puts "select atom1, id #{a.atom1}"
	tf.puts "select atom2, id #{a.atom2}"
	tf.puts "show spheres, atom1, atom2"
	tf.puts "color blue, atom1, atom2"
	tf.puts "draw_links atom1, atom2, red, red, #{Cylinder_diameter}, link_#{a.atom1}_#{a.atom2}"
end
	tf.puts 'zoom'
	tf.puts "png /home/ryanmt/Dropbox/coding/draw_er/#{file}.png"
	tf.puts "save /home/ryanmt/Dropbox/coding/draw_er/#{file}.pse,format=pse"
#	tf.puts 'log_close'
	tf.puts 'quit' 
	tf.flush
	%x[pymol #{tf.path}]
end




