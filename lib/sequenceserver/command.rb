# Methods of general use for command line applications
module SequenceServer
  class Command
    # Execute a command and return its stdout, stderr, and exit status.
    def self.execute(command)
      rfile = Tempfile.new('sequenceserver_result')
      efile = Tempfile.new('sequenceserver_error')
      [rfile, efile].each {|file| file.close}

      command = "#{command} > #{rfile.path} 2> #{efile.path}"
      system(command)
      status = $?.exitstatus

      return File.readlines(rfile.path), File.readlines(efile.path), status
    ensure
      rfile.unlink
      efile.unlink
    end
  end
end