#!/usr/bin/env ruby 
# test_ssequencehelpers.rb

require 'sequenceserver'
require 'lib/sequencehelpers'
require 'test/unit'



class Tester < Test::Unit::TestCase
  include SequenceServer::SequenceHelpers
  def test_guess_sequence_type_nucleotide
  #must 'correctly detect nucleotide sequence, even when it includes crap' do   
    ['AAAAAAAAAAAAAAAAAAAAAT',
     '           CAGATGCRRCAAAGCAAACGGCAA 34523453 652352',
     'ACCNNNNNNXXXXCAUUUUUU',
     "ACGT\n\t\t\nACCACGGACCACGAAAGCG"               ].each do |seq|
      assert_equal(:nucleotide, guess_sequence_type(seq), message="for #{seq}")
    end
  end

  def test_guess_sequence_type_aminoacid
  #must 'correctly detect aminoacid sequence, even when it includes a lot of crap' do
  ['ADSACGHKSJLFCVMGTL',
   '  345   KSSYPHYSPPPPHS      345 23453 652352',
   'GEYSNLNNNNNNXXXXSSSSSSSSSSSSSSSSSSSSSSS',
   "EE\n\t\t\n         \t\t\EEQRRQQSARTSRRQR"     ].each do |seq|
      assert_equal(:protein, guess_sequence_type(seq) , message="for #{seq}")
    end
  end

  def test_guess_sequence_type_impossible
    assert_equal(nil, guess_sequence_type('ACSFGT'), message='too little sequence')
  end

  ## Tests for type_of_sequences  (multi-fasta  kind of thing the user would enter)
  def test_type_of_sequences
    aa_multifasta = ">SDFDSF\nACGTGSDLKJGNLDKSJFLSDKJFLSDKOIU\n>asdfas\nasfasdfffffffffffffffffffff\n>alksjflkasdj slakdjf\nasdfasdfasdfljaslkdjf"
    aa_multifasta_including_short_seq_missing_lead = "ACGTGSDLKJGNLDKSJFLSDKJFLSDKOIU\n>asdfas\nasf\n>alksjflkasdj slakdjf\nasdfasdfasdfljaslkdjf"
    aa_singlesequence =  "ACGTGSDLKJGNLDKSJFLSDKJFLSDKOIU\n"
    nt_multifasta = ">asdf\nAAAAAAAAAAAAAAAAAAAAT\n>sfaslkj\nCAGATGCRRCAAAGCAAACGGCAA\n>asssssjlkj\nACCCANNNNNNXXXXCAUUUUUU"
    aa_nt_mix = ">alksjflkasdj slakdjf\nasdfasdfasdfljaslkdjf\n>ffffffassdf\nACGCNAGTGCCCCCCCCGANATGGGTGGTTXXXXXGGTG"

    assert_equal(:protein, type_of_sequences(aa_multifasta), 'aa_multifasta')
    assert_equal(:protein, type_of_sequences(aa_multifasta_including_short_seq_missing_lead ), 'aa_multifasta_short_seq_and_no_>')
    assert_equal(:protein, type_of_sequences(aa_singlesequence), 'single AA sequence')
    assert_equal(:nucleotide, type_of_sequences(nt_multifasta), 'nt_multifasta')
    assert_raise(ArgumentError, 'mixed aa and nt should raise') { type_of_sequences(aa_nt_mix) }
  end

  def test_sequence_type_to_blast_methods
    assert_equal ['blastp', 'tblastn'], blast_methods_for(:protein), 'blasts_for_protein'
    assert_equal ['blastn','tblastx','blastx'], blast_methods_for(:nucleotide), 'blasts_for_nucleotide'
    assert_equal ['blastp', 'tblastn','blastn','tblastx','blastx'], blast_methods_for(nil), 'blasts_for_nil'
  end

  def test_composition
    expected_comp = {"a"=>2, "d"=>3, "f"=>7, "s"=>3, "A"=>1}
    assert_equal(expected_comp, composition('asdfasdfffffAsdf'))
  end
  
  def test_construct_standard_sequence_hyperlink
    assert_equal "/get_sequence/:one/:abc def", construct_standard_sequence_hyperlink({:sequence_id => 'one', :databases => %w(abc def)})
    assert_equal nil, construct_standard_sequence_hyperlink({:sequence_id => ' one', :databases =>  %w(abc def)})
    assert_equal "/get_sequence/:MAL13P1.218/:abc def", construct_standard_sequence_hyperlink({:sequence_id => 'lcl|MAL13P1.218', :databases =>  %w(abc def)})
  end
end

class SystemHelperTester < Test::Unit::TestCase
  include SequenceServer::Helpers::SystemHelpers
  def test_process_advanced_blast_options
    assert_nothing_raised {process_advanced_blast_options('')}
    assert_nothing_raised {process_advanced_blast_options('-word_size 5')}
    assert_raise(ArgumentError, 'security advanced option parser'){process_advanced_blast_options('-word_size 5; rm -rf /')}
    assert_raise(ArgumentError, 'conflicting advanced option'){process_advanced_blast_options('-db roar')}
  end
end

class BlastTester < Test::Unit::TestCase
  include SequenceServer
  
  # Where the test data is stored. Might be changed in future or
  # made more available to other test classes
  def protein_db_dir
    "#{File.join(File.dirname(__FILE__), 'database', 'protein')}"
  end
  
  def setup
    # Create blast databases if they don't already exist
    Dir.chdir(protein_db_dir) do
      unless File.exist?('Sinvicta2-2-3.prot.subset.fasta.psi')
        `makeblastdb -in 'Sinvicta2-2-3.prot.subset.fasta' -parse_seqids`
      end
    end
  end
  
  def test_convert_blast_archive_to_hit_objects
    # qseqid qlen qstart qend evalue sseqid slen
    # PF14_0448 272 3 35  0.32  SI2.2.0_03512 1122
    # PF14_0448 272 223 252 1.6 SI2.2.0_02890 69
    # PF14_0448 272 225 253 2.8 SI2.2.0_08297 67
    new_hit_object = lambda do |qseqid, qlen, qstart, qend, evalue, sseqid, slen|
      h1 = Blast::Hit.new
      h1.qseqid = qseqid
      h1.qlen = qlen
      h1.qstart = qstart
      h1.qend = qend
      h1.evalue = evalue
      h1.sseqid = sseqid
      h1.slen = slen
      h1
    end
    
    hits = [
      new_hit_object.call('PF14_0448',272,3,35,0.32,'SI2.2.0_03512',1122),
      new_hit_object.call('PF14_0448',272,223,252,1.6,'SI2.2.0_02890',69),
      new_hit_object.call('PF14_0448',272,225,253,2.8,'SI2.2.0_08297',67),
    ]
    query = <<END
>PF14_0448 40S ribosomal protein S2, putative
MEDRGGFSRGFGRGVRGTRGRGGRGARGRGRGSAEDDLKNWVPVTKLGRLVKEGKIVSIEEIYLHSLPIKEYQIIDYFFQ
PNECSHPLKDDVVKIMPVQKQTRAGQRTRFKAFVAIGDGNGHCGLGVKCAKEVATAIRGAIISAKLSLIPVRRGYWGNKI
GDPHTVPMKVSGKCGSVRIRLVPAPRGTQIVGAPTTKKMLNFAGIKDCFSSSCGKTKTKGNFLRAIFNALSKTYGYLTPD
LWKVTNFDKSPYEEWSDFLETYQNLKGIKGTV
END
    b = Blast.blast_string_to_blast_archive('blastp', File.join(protein_db_dir, 'Sinvicta2-2-3.prot.subset.fasta'), query)
    actual_hits = b.convert_blast_archive_to_hit_objects('blast_formatter')
    assert_equal hits.length, actual_hits.length, "Length of array difference - #{hits} expected, #{actual_hits} found"
    hits.each_with_index do |hit, i|
      assert_equal hit, actual_hits[i]
    end
    
    q2 = <<END
MHPTVLDATGYTLLLSNFITFAILLWKGYKRKRKCPFYVLSLSDIFSATLMAVVLLVNHI
EAGIRLNYNWQNNTGGDMPNHTWTIQDKRFPFLQMHLREVDDLDVTLTCGMKDIFMHYGM
LLAALANAFTSLLTFAVQCNFNAAAIKRRCANVMKSSLKNAQLELPTDAEVKCLSELKST
SRERRIDAKVKSNVTQRIEKXXXXXXXXXXXXXXXXXXXXXXXXQLFVLSLAIYFLPILLSSILQMRGKHMCK
NTLAILRAKTNFTFTDGKKSQSRDSVEFTVPQGSRNDRSKTIIDVMKEGSCKENESIALE
IDRMVRTLDTIKLSLILCVLLWSPVFLGTLLRVYSCTRAPQWLTDVTFLSAILFGIVRNV
LNVNIVRIQEACTDANAKDNRIQPV
END
    hits = [
      new_hit_object.call('unnamed',398,1,200,6e-121,'SI2.2.0_09373',745),
      new_hit_object.call('unnamed',398,225,398,3e-101,'SI2.2.0_09373',745),
      new_hit_object.call('unnamed',398,260,310,0.008,'SI2.2.0_06144',491),
    ]
    b = Blast.blast_string_to_blast_archive('blastp', File.join(protein_db_dir, 'Sinvicta2-2-3.prot.subset.fasta'), q2)
    actual_hits = b.convert_blast_archive_to_hit_objects('blast_formatter')
    assert_equal hits.length, actual_hits.length, "Length of array difference - #{hits} expected, #{actual_hits} found"
    hits.each_with_index do |hit, i|
      assert hit==actual_hits[i]
    end
  end
end

class GraphicsTester < Test::Unit::TestCase
  include SequenceServer
  
  def test_evalue_to_lane_size
    g = Graphic.new
    assert_equal 2, g.evalue_to_lane_size(1)
    assert_equal 2, g.evalue_to_lane_size(0.1)
    assert_equal 6, g.evalue_to_lane_size(0.001)
    assert_equal g.evalue_to_lane_size(1e-50), g.evalue_to_lane_size(1e-51)
  end
end