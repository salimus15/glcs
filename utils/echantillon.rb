#!/usr/bin/env ruby

if ARGV.length != 1
	puts "Usage for all clean files #{ARGV[0]}: echantillon value "
	return	
end

files=Dir.glob("*.clean")

files.each do |file|
	puts file+"  ->  "+file.split('.').first+".ech"
	fin=File.open(file,'r')
	fout=File.open(file.split('.').first+".ech",'w')

	i=0
	fin.readlines.each do |line|
		if ((i).to_f.modulo(ARGV[0].to_f)).to_i==0
			fout.puts line
		end
		i=i+1
	end
	fin.close
	fout.close

end




