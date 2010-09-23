#!/usr/bin/ruby

require 'rubygems'
# require 'bio'
require 'yaml'
# require 'narray'

# Local (User) library
require 'ryan_lib'
require 'binary_search'

# Usage guidelines
if ARGV.size == 0
	puts "usage: #{File.basename(__FILE__)} < file>.yml ... "
	puts "output: <file>.yml ... "
	exit
end
total=[]
ARGV.each do |file|
	total <<  YAML::load(IO.read(file))
end
dist = []
dist_inter =[]
dist_intra =[]
inside = []
count_inter = 0
count_total = 0
intra_50=[]; inter_50=[];
total.compact!
total.delete(false)
total.each do |a|
	next if a.class ==! 'Array'
	a.each do |b|
#		next if b.class !='Hash'
		count_total += 1
		dist << b.fetch("dist")
		inside << b.fetch("inside")
		if b.fetch("inter")
			dist_inter << b.fetch("dist")
			count_inter += 1
		else
			dist_intra << b.fetch("dist")
		end
		if inside.last < 0.50
			intra_50 << dist_intra.last
			inter_50 << dist_inter.last
		end
	end
end
require 'simpler'
hist_inter = Simpler.new
hist_intra = Simpler.new
hist_global = Simpler.new
pie_inter = Simpler.new
double_hist = Simpler.new
range = (3..200).to_a
inter_perc = count_inter/count_total
intra_perc = (count_inter-count_total)/count_total
# Data:
=begin
dist_inter
dist_intra
dist
intra_perc
inter_perc	
=end
hist_inter.with(dist_inter) do |x|
	%Q{pdf(file="/home/ryanmt/Dropbox/coding/db_linker/inter.pdf", height=8, width=5)
	hist(#{x}, breaks=100,
	ylab = "Counts",
	xlab = "Length in Angstroms",
	main = "Histogram of Interchain Cross-Links"	
	)}
end
hist_intra.with(dist_intra) do |x|
	%Q{pdf(file="/home/ryanmt/Dropbox/coding/db_linker/intra.pdf", height=8, width=5)
	hist(#{x}, breaks=100,
	ylab = "Counts",
	xlab = "Length in Angstroms",
	main = "Histogram of Intrachain Cross-Links"
	)}
end
hist_global.with(dist) do |x|
	%Q{pdf(file="/home/ryanmt/Dropbox/coding/db_linker/global.pdf", height=5, width=5)
	hist(#{x}, breaks=100,
	ylab = "Counts",
	xlab = "Length in Angstroms",
	main = "Global Histogram of Cross-Links"
	)}
end
=begin hist_35.with(out_35) do |x|
	%Q{pdf(file="/home/ryanmt/Dropbox/coding/db_linker/h_35.pdf", height=5, width=5)
	hist(#{x}, breaks=100,
	ylab = "Counts",
	xlab = "Length in Angstroms",
	main = "Histogram of Cross-Links with < 35% insideness"
	)}
end
double_hist.with(intra_50.compact, inter_50.compact) do |intra, inter|
	%Q{require(plotrix)
	pdf(file="/home/ryanmt/Dropbox/coding/db_linker/double.pdf", height=5, width=5)
	input <- list(#{intra},#{inter})
	multihist(l, main="Histogram of Cross-Links with >= 50% outsideness", col=list("blue","green"))
	

} end
=end
# --------------------------
# Use the outputs and then plot them on a bar plot or a line plot... I don't see any other way to acheive this quickly!!!!
# --------------------------
double_hist.with(intra_50.compact, inter_50.compact) do |intra, inter|
	%Q{superhist <- function(x, filename= "/home/ryanmt/Dropbox/coding/db_linker/super_histograms.pdf",
dev = "pdf", title = "Red for Intrachain X-links, Blue for Interchain X-links, of those with >= 50% outsideness", nbreaks ="Sturges") {
junk = NULL
grouping = NULL
for(i in 1:length(x)) {
junk = c(junk,x[[i]])
grouping <- c(grouping, rep(i,length(x[[i]]))) }
grouping <- factor(grouping)
n.gr <- length(table(grouping))
xr <- range(junk)
histL <- tapply(junk, grouping, hist, breaks=nbreaks, plot = FALSE)
maxC <- max(sapply(lapply(histL, "[[", "counts"), max))
if(dev == "pdf") { pdf(filename, version = "1.4") } else{}
if((TC <- transparent.cols <- .Device %in% c("pdf", "png"))) {
cols <- hcl(h = seq(30, by=360 / n.gr, length = n.gr), l = 65, alpha = 0.5) }
else {
h.den <- c(10, 15, 20)
h.ang <- c(45, 15, -30) }
if(TC) {
plot(histL[[1]], xlim = xr, ylim= c(0, maxC), col = cols[1], xlab = "Distance in Angstroms", main = title) }
else { plot(histL[[1]], xlim = xr, ylim= c(0, maxC), density = h.den[1], angle = h.ang[1], xlab = "Distance in Angstroms") }
if(!transparent.cols) {
for(j in 2:n.gr) plot(histL[[j]], add = TRUE, density = h.den[j], angle = h.ang[j]) } else {
for(j in 2:n.gr) plot(histL[[j]], add = TRUE, col = cols[j]) }
invisible()
if( dev == "pdf") {
dev.off() }
}
	input=list(#{intra},#{inter})
d1 = rnorm(1:100)
d2 = rnorm(1:100) + 4
# the input object MUST be a list!
l1 = list(d1,d2)
superhist(input, nbreaks="Sturges")	}
# Red , Blue
	#ylab = "Counts",
	#xlab = "Length in Angstroms",
	#main = "Histogram of Cross-Links with < 65% insideness"
end
count_intra=(count_total-count_inter)
bar_in = [count_inter, count_intra]
pie_inter.with(bar_in) do |bar|
	%Q{pdf(file="/home/ryanmt/Dropbox/coding/db_linker/pied.pdf", height=8, width=5)
names(#{bar}) <- c("Interchain X-links", "Intrachain X-links")
	barplot(#{bar}, col=c("purple", "green3"))}
end
hist_inter.show!
hist_intra.show!
hist_global.show!
pie_inter.show!
double_hist.show!
# Linkers output to yaml	
	ymlfile1 = "out.yml"
	puts "creating: #{ymlfile1}" if $VERBOSE
	File.open(ymlfile1, "w") {|out| out.print puts "dist: #{dist.size}\n inside: #{inside.size}\n inter/total: #{count_inter}/#{count_total} or #{count_inter/count_total.to_f}\n inter distances: #{dist_inter.size}\n intra distances: #{dist_intra.size}" }
# Inside added to database file
	#~ out2file = "outs.yml",
	#~ short_out = linkers.map do |a| a.to_archive end 
	#~ File.open(out2file, "w") {|output| output.print short_out.to_yaml }
