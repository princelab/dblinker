#!/usr/bin/ruby

require 'rubygems'
# require 'bio'
require 'yaml'
# require 'narray'

# Local (User) library
# require 'ryan_lib'
# require 'binary_search'

# Usage guidelines
if ARGV.size == 0
	puts "usage: #{File.basename(__FILE__)} < file>.yml ... "
	puts "output: <file>.yml ... "
	exit
end
total=[]
ARGV.each do |file|
	puts 'file: ' + file.to_s
	total <<  YAML::load(IO.read(file))
end
dist = []
total.compact!
total.delete(false)
total.each do |a|
	next if a.class ==! 'Array'
	a.each do |b|
		dist << b.fetch("dist_to_solvent_1")
		dist << b.fetch("dist_to_solvent_2")
	end
end
require 'simpler'
hist_dist = Simpler.new
# range = (3..200).to_a
hist_dist.with(dist) do |x|
	%Q{pdf(file="/home/ryanmt/Dropbox/coding/dblinker/spec/dist.pdf", height=8, width=5)
	hist(#{x}, breaks=100,
	ylab = "Counts",
	xlab = "Length in Angstroms",
	main = "Histogram of solvent accessibility distances for all 'NZ' atoms"
	)}
end
hist_dist.show!
