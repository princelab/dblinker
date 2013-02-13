#!/usr/bin/env ruby
# require 'putsv' # a gem to provide my traditional tools... Great idea Ryan!

def putsv(str)
  puts str if $VERBOSE
end

require 'yaml'
Cut_off=0.8
Cylinder_diameter=0.1
PYMOL_RAY_DPI = 1200
PYMOL_DPI = 300
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
  dirname = File.absolute_path(File.dirname(file))
	filename = File.basename(file, ".xlink.yml")
  outfile = File.join(dirname, filename)
  %w{.pse .png}.map.with_index {|a,i| puts "Outfile ##{i+1}: #{outfile+a}"}
	%x[rm tmp.pml] if File.exists?('tmp.pml')
	tmp_file = File.open('tmp.pml','w+')
#	tmp_file.puts 'log_open ("tmp_log.txt")'
	tmp_file.puts 'run draw_links.py'
	tmp_file.puts "load #{filename}.pdb"
	tmp_file.puts 'hide everything'
	tmp_file.puts 'show cartoon'
outs.each do |a|
putsv a.atom2
putsv a.atom1
	tmp_file.puts "select atom1, id #{a.atom1}"
	tmp_file.puts "select atom2, id #{a.atom2}"
	tmp_file.puts "show spheres, atom1, atom2"
	tmp_file.puts "color blue, atom1, atom2"
	tmp_file.puts "draw_links atom1, atom2, red, red, #{Cylinder_diameter}, link_#{a.atom1}_#{a.atom2}"
end
	tmp_file.puts 'zoom'
  tmp_file.puts "ray #{PYMOL_RAY_DPI},#{PYMOL_RAY_DPI}"
	tmp_file.puts "png #{outfile}.png, dpi=#{PYMOL_DPI}"
	tmp_file.puts "save #{outfile}.pse,format=pse"
#	tmp_file.puts 'log_close'
	tmp_file.puts 'quit' 
	tmp_file.flush
  putsv %x[pymol -c #{tmp_file.path}]
end




