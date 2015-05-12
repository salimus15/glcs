#!/usr/bin/env ruby







class GraphMat
   def initialize(filename,fout)
	puts "Processing : "+filename
	@filename=filename
	@fout=fout
      	@name=filename.split('.').first
      	@param_l = @name.scan(/\d+/).at(4)
      	@param_k = @name.scan(/\d+/).at(3)
      	@param_latency = @name.scan(/\d+/).at(5)
   end

   def export(type)
	puts("Exporting figure"+" #{@name}_#{type}.png k = #{@param_k}, l = #{@param_l}, latency = #{@param_latency}")
	@fout.puts("load #{@filename}")
	@fout.puts("f = figure('Visible','off');")
	@fout.puts("figure('Name','#{@name.split('_').first} matrix','NumberTitle','off');")
	self.figure_nocolor(type)
	@fout.puts("legend('k = #{@param_k}, l = #{@param_l}, latency = #{@param_latency}');")
	@fout.puts("set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 50 50]);")
	@fout.puts("print('-dpng', '#{@name}_#{type}.png','-zbuffer');")
	@fout.puts("close(findobj('type','figure','Name','#{@name.split('_').first} matrix'))")
   end

   def figure_color(type,color)
	@fout.puts("load #{@filename}")
	@fout.puts("#{type}(#{@name},#{color}); hold on;")
   end
 
   def figure_nocolor(type)
	@fout.puts("load #{@filename}")
	@fout.puts("#{type}(#{@name}); hold on;")
   end
 

   def label()
	return "#{@name.split('_').first} : k = #{@param_k}, l = #{@param_l}, latency = #{@param_latency}"
   end

end




#read clean files
files=Dir.glob("*.ech")
fout=File.open("graphClean.m","w")

#process files
files.each do |file|
	g=GraphMat.new(file,fout)
	g.export("loglog")
	g.export("plot")
	g.export("semilogx")
	g.export("semilogy")
end
	
#prepare plot color gradient output
fout.puts("figure('Name','Comparison','NumberTitle','off');")
fout.puts("hold all")
fout.puts("couleurs = hsv(#{files.length})")
fout.puts("set(gca,'colororder',couleurs)")


labels=[]	
files.each do |file|
	g=GraphMat.new(file,fout)
	g.figure_nocolor("semilogx")
	labels.push(g.label)
end
	
fout.puts("set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 50 50]);")
fout.puts("print('-dpng', 'loglog.png','-zbuffer');")
fout.puts("close(findobj('type','figure','Name','Comparison'))")

#for i in 0..files.length-1
#	fout.puts("h(#{i+1})=semilogx(#{files[i].to_s.split('.').first});hold on;")
	#fout.puts("legend('#{file.split('.').first}');")
#end 
#fout.puts("legend(h,num2str([1:#{files.length}].'));")

fout.close


