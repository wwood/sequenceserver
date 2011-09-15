module SequenceServer
  class Graphic
    # Given a query_id and an Array of SequenceServer::Blast::Hit objects
    # create a graphic, and remove the hits from the array that pertain
    # to the hit.
    #
    # TODO for blast graphic2 branch 
    #* thinner gene sizes, particularly when there is many hits
    #* add hyperlinks to the specific query results
    #* make the different frames different colours for blastx, tblastn and tblastx
    #* check for the possible bug where there is no hits for the first sequence (hits won't be returned in the blast_formatter ryo response)
    #* check for bug where the query_id doesn't match the whole of the query name
    def format_query_graphical_overview(query_id, hits, query_length)
      # extract the hits for this query from the entire hits array
      pertinent_hits = []
      hits.each_with_index do |hit,i|
        next if hit.nil?
        if hit.qseqid == query_id
          pertinent_hits.push hit
          hits[i] = nil
        else
          break
        end
      end
      # Each canvas name should be unique
      canvas_name = "canvas_#{rand(100000).round}"
      # make the canvas height big enough to fit all the hits
      canvas_height = 50+20*pertinent_hits.length
      line = "<div id=\"container\"><canvas id=\"#{canvas_name}\" width='5000' height='#{canvas_height}'></canvas></div>\n"
      # Hit instance_variables = [:qseqid, :qlen, :qstart, :qend, :evalue, :sseqid, :slen]
      line += "<script>\nvar canvas = document.getElementById('#{canvas_name}');\n";
      line += "chart1 = new Scribl(canvas, 650);\n"
      # Force the canvas scale to not zoom in, as this isn't obvious enough for the casual user 
      line += "chart1.scale.min=1;\n"
      line += "chart1.scale.max=#{query_length};\n"
      line += "chart1.scale.auto = false;" #don't want intelligent start and stops
      line += "chart1.laneSizes = 18;\n" #use smaller than default track sizes by default.
      # Add Genes      position, length, orientation
      pertinent_hits.each_with_index do |hit, i|
        line += "gene#{i} = chart1.addFeature( new Rect('hit',#{hit.qstart},#{hit.qend-hit.qstart},'+'));\n"
      end
      line += "chart1.draw();\n</script>\n"
    end
  end
end