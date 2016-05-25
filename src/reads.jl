
function make_fqparser( filename; forcegzip=false )
   if isgzipped( filename ) || forcegzip
      fopen = open( filename, "r" ) 
      to_open = ZlibInflateInputStream( fopen, reset_on_end=true )
   else
      to_open = filename
   end 
   open( to_open, FASTQ )
end

function read_chunk!( chunk, parser )
   i = 1
   while i <= length(chunk) && read!( parser, chunk[i] )
      i += 1
   end
   while i <= length(chunk)
      pop!(chunk) # clean up if we are at the end
   end
   parser
end

function allocate_chunk( parser; size=100000 )
  chunk = Vector{eltype(parser)}( size )
  for i in 1:length(chunk)
     chunk[i] = eltype(parser)()
  end
  chunk
end

allocate_rref( size=250000; rreftype=RemoteRef{Channel{Any}} ) = Vector{rreftype}( size )

function resize_rref!( rref, subsize )
   @assert subsize <= length(rref)
   while subsize <= length(rref)
      pop!(rref) # clean up if we are at the end
   end
end

function sendto(p::Int, nm, val)
   ref = @spawnat(p, eval(Main, Expr(:(=), nm, val)))
end

macro sendto(p, nm, val)
   return :( sendto($p, $nm, $val) )
end

macro broadcast(nm, val)
   quote
      @sync for p in workers()
         @async sendto(p, $nm, $val)
      end
   end
end

function get_path( offset, leng, nodelist, sg )
   used = 0
   path = IntSet()
   len = leng
   has_started = false
   for n in nodelist
      n > length(sg.nodelen) && break
      if used <= offset < used+sg.nodelen[n]
         len -= (used+sg.nodelen[n]) - offset
         has_started = true
         push!(path, n)
      elseif has_started
         if len > 0
            push!(path, n)
            len -= sg.nodelen[n]
         else
            break
         end
      end
      used += sg.nodelen[n]
   end
   path
end

function cmp_paths( path::IntSet, aln::SGAlignment )

   partial = true
   complete = true
   if length(path) > 0 && (aln.path[1].node == first(path))
      shift!(path)
      for i in 2:length(aln.path)
         if length(path) == 0 || aln.path[i].node != first(path)
            complete = false
            break
         end
         shift!(path)
      end
   else
   #   println("path: $(best.path) vs. $path, $(read.name)")
      partial = false
      complete = false
   end
   partial,complete
end

# Use this function for true positive rates in the accuracy branch.
function ismappedcorrectly( read::SeqRecord, avec::Vector{SGAlignment}, lib::GraphLib )
   failed = false
   ord    = sortperm( avec, by=score )
   best   = avec[ ord[end] ]
   gene   = best.path[1].gene
   spl    = split( read.name, '/' )[end] |> x->split( x, ';' )
   @assert( length(spl) == 3, "ERROR: Incorrect format for simulated read name, $(spl)!" )
   nodes  = split( spl[1], '_' )[end] |> x->split( x, '-' )
   offset_m1 = split( spl[2], ':' )[end] |> x->split( x, '-' )
   offset_m2 = split( spl[3], ':' )[end] |> x->split( x, '-' )
   #@assert( length(offset) == 2, "ERROR: Incorrect format for simulated read name, $(read.name)!" )
   if length(offset_m1) != 2 || length(offset_m2) != 2
      return false,false,true
   end
   nlist = map( x->parse(Int, x), nodes )

   off_m1 = parse(Int, offset_m1[1]), parse(Int, offset_m1[2])
   mate1 = get_path( off_m1[1], off_m1[2] - off_m1[1], nlist, lib.graphs[gene] )
   part_m1,comp_m1 = cmp_paths( mate1, best )

   off_m2 = parse(Int, offset_m2[1]), parse(Int, offset_m2[2])
   mate2 = get_path( off_m2[1], off_m2[2] - off_m2[1], nlist, lib.graphs[gene] )
   part_m2,comp_m2 = cmp_paths( mate2, best )

   (part_m1 || part_m2),(comp_m1 || comp_m2),failed
end

process_reads!( parser, param::AlignParam, lib::GraphLib,
                quant::GraphLibQuant, multi::Vector{Multimap}; 
                bufsize=50, sam=false, simul=false) = _process_reads!( parser, param, lib, quant,
                                                      multi, bufsize=bufsize, sam=sam, simul=simul )

function _process_reads!( parser, param::AlignParam, lib::GraphLib, quant::GraphLibQuant, 
                         multi::Vector{Multimap}; bufsize=50, sam=false, simul=false )
  
   const reads  = allocate_chunk( parser, size=bufsize )
   mean_readlen = 0.0
   total        = 0
   mapped       = 0
   correct_part = 0
   correct_full = 0
   if sam
      stdbuf = BufferedOutputStream( STDOUT )
      write_sam_header( stdbuf, lib )
   end
   while length(reads) > 0
      read_chunk!( reads, parser )
      total += length(reads)
      for i in 1:length(reads)
         align = ungapped_align( param, lib, reads[i] )
         if !isnull( align )
            if length( align.value ) > 1
               push!( multi, Multimap( align.value ) )
            else
               count!( quant, align.value[1] )
               sam && write_sam( stdbuf, reads[i], align.value[1], lib )
            end
            if simul
               part,full,fail = ismappedcorrectly( reads[i], align.value, lib )
               correct_part += part ? 1 : 0
               correct_full += full ? 1 : 0
               !fail && (mapped += 1)
            else
               mapped += 1
            end
            @fastmath mean_readlen += (length(reads[i].seq) - mean_readlen) / mapped
         end
      end
   end # end while
   if sam
      close(stdbuf)
   end
   mapped,correct_part,correct_full,total,mean_readlen
end

# paired end version
process_paired_reads!( fwd_parser, rev_parser, param::AlignParam, lib::GraphLib,
                quant::GraphLibQuant, multi::Vector{Multimap}; 
                bufsize=50, sam=false) = _process_paired_reads!( fwd_parser, rev_parser, param, lib, quant,
                                                                 multi, bufsize=bufsize, sam=sam )

function _process_paired_reads!( fwd_parser, rev_parser, param::AlignParam, lib::GraphLib, quant::GraphLibQuant,
                                 multi::Vector{Multimap}; bufsize=50, sam=false )

   const fwd_reads  = allocate_chunk( fwd_parser, size=bufsize )
   const rev_reads  = allocate_chunk( rev_parser, size=bufsize )
   mean_readlen = 0.0
   total        = 0
   mapped       = 0
   if sam
      stdbuf = BufferedOutputStream( STDOUT )
      write_sam_header( stdbuf, lib )
   end
   while length(fwd_reads) > 0 && length(rev_reads) > 0
      read_chunk!( fwd_reads, fwd_parser )
      read_chunk!( rev_reads, rev_parser )
      total += length(fwd_reads)
      for i in 1:length(fwd_reads)
         fwd_aln,rev_aln = ungapped_align( param, lib, fwd_reads[i], rev_reads[i] )
         if !isnull( fwd_aln ) && !isnull( rev_aln )
            if length( fwd_aln.value ) > 1
               push!( multi, Multimap( fwd_aln.value ) )
               push!( multi, Multimap( rev_aln.value ) )
            else
               count!( quant, fwd_aln.value[1], rev_aln.value[1] )
               sam && write_sam( stdbuf, fwd_reads[i], fwd_aln.value[1], lib, 
                                 paired=true, first=true, is_pair_rc=param.is_pair_rc )
               sam && write_sam( stdbuf, rev_reads[i], rev_aln.value[1], lib, 
                                 paired=true, first=false, is_pair_rc=param.is_pair_rc )
            end
            mapped += 1
            @fastmath mean_readlen += (length(fwd_reads[i].seq) - mean_readlen) / mapped
         end
      end
   end # end while
   if sam
      close(stdbuf)
   end
   mapped,total,mean_readlen
end
