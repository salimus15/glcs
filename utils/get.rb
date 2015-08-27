#!/usr/bin/env ruby


result = []

files=Dir.glob("*.txt")
files.each do |file|
  puts "Processing "+file
  File.open(file,"r") do |handle|
    handle.each_line do |line|
	  if(line.include?('KSP preconditioned resid norm'))	
      result << line.split(' ').last
	  end
    end
  end # closes automatically when EOF reached
  result.flatten!
  
  
  puts file.split('.').first+".clean  records : "+result.length.to_s
  fout=File.open(file.split('.').first+".clean","w")
  fout.puts result
  fout.close
  
  result.clear
  
end

