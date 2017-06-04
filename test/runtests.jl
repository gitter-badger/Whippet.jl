using Base.Test

using DataStructures
using BufferedStreams
using Bio.Seq
using FMIndexes
using IntArrays
using IntervalTrees
using Libz
using Distributions
using Requests

include("../src/types.jl")
include("../src/timer.jl")
include("../src/sgkmer.jl")
include("../src/fmindex_patch.jl")
include("../src/refset.jl")
include("../src/graph.jl")
include("../src/edges.jl")
include("../src/index.jl")
include("../src/align.jl")
include("../src/quant.jl")
include("../src/reads.jl")
include("../src/ebi.jl")
include("../src/paired.jl")
include("../src/events.jl")
include("../src/io.jl")
include("../src/diff.jl")

@testset "SG Kmers & Seq" begin
   @test sgkmer(dna"ATG") == Bio.Seq.DNAKmer(dna"ATG")
   @test sgkmer("ATG") == Bio.Seq.DNAKmer(dna"ATG")
   @test isa(sgkmer(dna"ATG"), SGKmer{3})
   @test kmer_index( sgkmer(dna"ATG") ) == kmer_index( Bio.Seq.DNAKmer(dna"ATG") )
   @test Int(kmer_index_trailing( dna"ATG" ))+1 == kmer_index( Bio.Seq.DNAKmer(dna"ATG") )
   @test Int(kmer_index_trailing( dna"ATGR" )) == 0
   seq = dna""
   for i in 1:32
      seq = i % 2 == 0 ? seq * dna"C" : seq * dna"T"
      curkmer = sgkmer(seq)
      @test isa( curkmer, SGKmer{i} )
   end
end
@testset "Splice Graphs" begin
   gtf = IOBuffer("# gtf file test
chr0\tTEST\texon\t6\t20\t.\t+\t.\tgene_id \"one\"; transcript_id \"def\";
chr0\tTEST\texon\t31\t40\t.\t+\t.\tgene_id \"one\"; transcript_id \"def\";
chr0\tTEST\texon\t54\t62\t.\t+\t.\tgene_id \"one\"; transcript_id \"def\";
chr0\tTEST\texon\t76\t85\t.\t+\t.\tgene_id \"one\"; transcript_id \"def\";
chr0\tTEST\texon\t6\t40\t.\t+\t.\tgene_id \"one\"; transcript_id \"int1_alt3\";
chr0\tTEST\texon\t51\t62\t.\t+\t.\tgene_id \"one\"; transcript_id \"int1_alt3\";
chr0\tTEST\texon\t76\t85\t.\t+\t.\tgene_id \"one\"; transcript_id \"int1_alt3\";
chr0\tTEST\texon\t6\t20\t.\t+\t.\tgene_id \"one\"; transcript_id \"apa_alt5\";
chr0\tTEST\texon\t31\t40\t.\t+\t.\tgene_id \"one\"; transcript_id \"apa_alt5\";
chr0\tTEST\texon\t54\t65\t.\t+\t.\tgene_id \"one\"; transcript_id \"apa_alt5\";
chr0\tTEST\texon\t76\t90\t.\t+\t.\tgene_id \"one\"; transcript_id \"apa_alt5\";
chr0\tTEST\texon\t11\t20\t.\t+\t.\tgene_id \"single\"; transcript_id \"ex1_single\";
chr0\tTEST\texon\t11\t20\t.\t-\t.\tgene_id \"single_rev\"; transcript_id \"single_rev\";
chr0\tTEST\texon\t11\t20\t.\t-\t.\tgene_id \"kissing\"; transcript_id \"def_kiss\";
chr0\tTEST\texon\t21\t30\t.\t-\t.\tgene_id \"kissing\"; transcript_id \"def_kiss\";
chr0\tTEST\texon\t11\t30\t.\t-\t.\tgene_id \"kissing\"; transcript_id \"ret_kiss\";
")

   flat = IOBuffer("# refflat file test (gtfToGenePred -genePredExt test.gtf test.flat)
def\tchr0\t+\t5\t85\t85\t85\t4\t5,30,53,75,\t20,40,62,85,\t0\tone\tnone\tnone\t-1,-1,-1,-1,
int1_alt3\tchr0\t+\t5\t85\t85\t85\t3\t5,50,75,\t40,62,85,\t0\tone\tnone\tnone\t-1,-1,-1,
apa_alt5\tchr0\t+\t5\t90\t90\t90\t4\t5,30,53,75,\t20,40,65,90,\t0\tone\tnone\tnone\t-1,-1,-1,-1,
ex1_single\tchr0\t+\t10\t20\t10\t20\t1\t10,\t20,\t0\tsingle\tnone\tnone\t-1,
")

   gtfref  = load_gtf( gtf )
   flatref = load_refflat( flat )

   @testset "Gene Annotation" begin

      for gene in keys(flatref)
         @test gtfref[gene].don    == flatref[gene].don
         @test gtfref[gene].acc    == flatref[gene].acc
         @test gtfref[gene].txst   == flatref[gene].txst
         @test gtfref[gene].txen   == flatref[gene].txen
         @test gtfref[gene].length == flatref[gene].length
         for i in 1:length(flatref[gene].reftx)
            flattx = flatref[gene].reftx[i]
            gtftx  = gtfref[gene].reftx[i]
            @test flattx.don == gtftx.don
            @test flattx.acc == gtftx.acc
         end
      end

   end

                               # fwd     rev
   buffer1   = dna"AAAAA"      # 1-5     96-100
   utr5      = dna"TTATT"      # 6-10    91-95
   exon1     = dna"GCGGATTACA" # 11-20   81-90
   int1      = dna"TTTTTTTTTT" # 21-30   71-80
   exon2     = dna"GCATTAGAAG" # 31-40   61-70
   int2      = dna"GGGGGGGGGG" # 41-50   51-60
   exon3alt3 = dna"CCT"        # 51-53   48-50
   exon3def  = dna"CTATGCTAG"  # 54-62   39-47
   exon3alt5 = dna"TTC"        # 63-65   36-38
   int3      = dna"CCCCCCCCCC" # 66-75   26-35
   exon4     = dna"TTAGACAAGA" # 76-85   16-25
   apa       = dna"AATAA"      # 86-90   11-15
   buffer2   = dna"AAAAAAAAAA" # 91-100  1-10

   fwd = buffer1 * utr5 * exon1 * int1 * 
         exon2 * int2 * 
         exon3alt3 * exon3def * exon3alt5 * int3 * 
         exon4 * apa * buffer2
   rev = reverse_complement(fwd)

   genome = fwd * rev

   expected_one = utr5 * exon1 * 
               int1 *
               exon2 *
               exon3alt3 * 
               exon3def *  
               exon3alt5 *
               exon4 * 
               apa 

   expected_sin = dna"GCGGATTACA"
   expected_kis = dna"AAAAAAAAAATGTAATCCGC"

   kmer_size = 2 # good test size

   graph_one = SpliceGraph( gtfref["one"], genome, kmer_size )
   graph_sin = SpliceGraph( gtfref["single"], genome, kmer_size )
   graph_rev = SpliceGraph( gtfref["single_rev"], genome, kmer_size )
   graph_kis = SpliceGraph( gtfref["kissing"], genome, kmer_size )

   @testset "Graph Building" begin
      @test graph_one.seq == expected_one
      @test graph_sin.seq == expected_sin
      @test graph_kis.seq == expected_kis

      @test length(graph_one.annopath) == length(gtfref["one"].reftx)
      @test graph_one.annopath[1] == IntSet([1,3,5,7]) # path of def
      @test graph_one.annopath[2] == IntSet([1,2,3,4,5,7]) # path of int1_alt3
      @test graph_one.annopath[3] == IntSet([1,3,5,6,7,8]) # path of apa_alt5
   end

   # Build Index (from index.jl)
   xcript  = dna""
   xoffset = Vector{UInt64}()
   xgenes  = Vector{GeneName}()
   xinfo   = Vector{GeneInfo}()
   xlength = Vector{Float64}()
   xgraph  = Vector{SpliceGraph}()

   runoffset = 0

   for g in keys(gtfref)
      println(STDERR, g)
      curgraph = SpliceGraph( gtfref[g], genome, kmer_size )
      xcript  *= curgraph.seq
      push!(xgraph, curgraph)
      push!(xgenes, g)
      push!(xinfo, gtfref[g].info )
      push!(xlength, gtfref[g].length )
      push!(xoffset, runoffset)
      runoffset += length(curgraph.seq)   
   end

   fm = FMIndex(twobit_enc(xcript), 4, r=1, program=:SuffixArrays, mmap=true)

   edges = build_edges( xgraph, kmer_size )

   #println(edges.left)
   #println(edges.right)

   lib = GraphLib( xoffset, xgenes, xinfo, xlength, xgraph, edges, fm, true, kmer_size )

   @testset "Kmer Edges" begin
      left  = [dna"CA", dna"AG", dna"AG", dna"TC", dna"AA"]
      right = [dna"GC", dna"CC", dna"CT", dna"TT", dna"TG"]
      lkmer = map( x->kmer_index(SGKmer(x)), left )
      rkmer = map( x->kmer_index(SGKmer(x)), right )
      #println(lkmer[1])
      #println(rkmer[1])
      for i in 1:4^kmer_size
         if i in lkmer
            @test isdefined(edges.left, i)
            @test typeof(edges.left[i]) == Vector{SGNode}
            @test issorted(edges.left[i], lt=sortlt)
         else
            @test !isdefined(edges.left, i)
         end
         if i in rkmer
            @test isdefined(edges.right, i)
            @test typeof(edges.right[i]) == Vector{SGNode}
            @test issorted(edges.right[i], lt=sortlt)
         else
            @test !isdefined(edges.right, i)
         end
      end
      exon1_lind = lkmer[1]
      exon2_rind = rkmer[1]
      @test intersect( edges.left[exon1_lind], edges.right[exon2_rind] ) == edges.right[exon2_rind]
      @test edges.left[exon1_lind] ∩ edges.right[exon2_rind] == edges.right[exon2_rind]
   end

   @testset "Saving and Loading Index" begin
      println(STDERR, "Saving test index...")
      println(STDERR, lib)
      @timer open("test_index.jls", "w+") do io
         serialize( io, lib )
      end
      println(STDERR, "Loading test index...")
      @timer lib = open(deserialize, "test_index.jls")
      println(STDERR, lib)
   end

   @testset "Alignment and SAM Format" begin
      # reads
      fastq = IOBuffer("@1S10M%11,20%exon1
NGCGGATTACA
+
#BBBBBBBBBB
@1S9M%54,62%exon3def
NCTATGCTAG
+
#BBBBBBBBB
@1S15M%51,65%alt3-exon3-alt5
NCCTCTATGCTAGTTC
+
#BBBBBBBBBBBBBBB
@11M10N10M%10,40%exon1-exon2
NGCGGATTACAGCATTAGAAG
+
#BBBBBBBBBBBBBBBBBBBB
@5M10N6M%16,36%exon1trunc-exon2trunc
TTACAGCATTN
+
BBBBBBBBBB#
@10M33N9M1S%11,62%exon1-exon3def
GCGGATTACACTATGCTAGN
+
BBBBBBBBBBBBBBBBBBB#
@10M33N9M%11,62%exon1-exon3def:rc
CTAGCATAGTGTAATCCGC
+
BBBBBBBBBBBBBBBBBBB
@10M55N10M1S%11,85%exon1-exon4full
GCGGATTACATTAGACAAGAN
+
BBBBBBBBBBBBBBBBBBBB#
@11M55N4M%10,78%exon1-exon4_4bp
NGCGGATTACATTAG
+
#IIIIIIIIIIIIII
@2M55N10M%19,85%exon1_2bp-exon4:rc
TCTTGTCTAATG
+
IIIIIIIIIIII
")

      score_range = 0.05
      param = AlignParam( 1, 2, 4, 4, 4, 5, 1, 2, 1000, score_range, 0.7,
                          false, false, true, false, true )
      quant = GraphLibQuant( lib )
      multi = Vector{Multimap}()

      typealias DNASeqType Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{2}}
      fqparse = FASTQReader{DNASeqType}( BufferedInputStream(fastq), Bio.Seq.ILLUMINA18_QUAL_ENCODING, DNA_C )
      reads  = allocate_chunk( fqparse, size=10 )
      read_chunk!( reads, fqparse )

      @test length(reads) == 10
      for r in reads
         println(STDERR, r)

         align = ungapped_align( param, lib, r )
         println(STDERR, align)

         flush(STDERR)
         @test !isnull( align )
         @test length( align.value ) >= 1
         @test all(map( x->x.isvalid, align.value))
         scores = map( x->identity(x, length(r.seq)), align.value )
         @test maximum(scores) - minimum(scores) <= score_range

         if length(search(r.name, ":rc")) > 0
            @test align.value[1].strand == false
         else
            @test align.value[1].strand == true
         end

         best_ind = indmax(scores)

         ex_num = length(split(r.name, '-', keep=false))
         @test length(align.value[best_ind].path) == ex_num

         count!( quant, align.value[best_ind] )

         const curgene   = align.value[best_ind].path[1].gene
         const firstnode = align.value[best_ind].path[1].node
         const lastnode  = align.value[best_ind].path[end].node
         const curgraph  = lib.graphs[ curgene ]

         cigar,positions = split(r.name, '%', keep=false)[1:2]
         first,last   = split(positions, ',', keep=false) |> y->map(x->parse(Int,x), y)
         # Test SAM Format
         @test cigar == cigar_string( align.value[best_ind], curgraph, true, length(r.seq) )[1]
         test_cigar,endpos = cigar_string( align.value[best_ind], curgraph, true, length(r.seq) )
         #println(STDERR, "cigar = $test_cigar")
         @test first == Int(curgraph.nodecoord[firstnode] + (align.value[best_ind].offset - curgraph.nodeoffset[firstnode]))
         #@test last  == Int(curgraph.nodecoord[lastnode] + (endpos - curgraph.nodeoffset[lastnode]))
         # test readlength = number of M and S entries in cigar  
         # test SAM offset is correct for both '+' and '-' genes.
         # test that cigar reversal works for '-' strand genes.
      end 

      @testset "Quantification" begin
         calculate_tpm!( quant, readlen=20 )

         #println(quant)
      end

      @testset "Event Building" begin
         out = IOBuffer(true,true)
         bs  = BufferedOutputStream(out)
         output_psi_header( bs )
         for g in 1:length(lib.graphs)
            name = lib.names[g]
            chr  = lib.info[g].name
            strand = lib.info[g].strand ? '+' : '-'
            _process_events( bs, lib.graphs[g], quant.quant[g], (name,chr,strand), isnodeok=false )
         end
         flush(bs)
         seek(out,0)
         for l in eachline(out)
            spl = split( l, '\t' )
            @test length(spl) == 13
            print(STDERR,l)
         end
      end

      @testset "EBI Accessions & HTTP Streaming" begin
         ebi_res = ident_to_fastq_url("SRR1199010") # small single cell file
         @test ebi_res.success
         @test !ebi_res.paired

         println(STDERR, "Streaming fastq file from $(ebi_res.fastq_1_url)")
         parser, response = make_http_fqparser( "http://" * ebi_res.fastq_1_url )
         reads  = allocate_chunk( parser, size=50 )
         cnt    = 0
         while length(reads) > 0
            read_http_chunk!( reads, parser, response )
            cnt += length(reads)
         end
         @test cnt == 128482 # correct number of reads in file
         @test typeof(run(`julia ../bin/whippet-quant.jl --ebi SRR1199010 -x test_index`)) <: Void
      end
   end
end
