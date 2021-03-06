#=
Basic Types and Aliases
=#

const Mb = 1_000_000
const GENOMESIZE = 3235Mb

# ALIAS NEW TYPES FOR INCREASED CODE READABILITY
if GENOMESIZE < typemax(UInt32)
   typealias CoordInt UInt32
else
   typealias CoordInt UInt64
end

typealias CoordTuple Tuple{Vararg{CoordInt}}
typealias CoordArray Vector{CoordInt}
typealias CoordTree  IntervalTree{CoordInt,IntervalTrees.Interval{CoordInt}}
typealias CoordSet   Set{Tuple{CoordInt,CoordInt}}

typealias ExonInt    UInt16
typealias ExonTree   IntervalTree{ExonInt,IntervalTrees.Interval{ExonInt}}

typealias SGSequence Bio.Seq.BioSequence{DNAAlphabet{4}}

typealias GeneName    String
typealias SeqName     String
typealias RefSeqId    String
typealias ASCIIString String

typealias GeneMeta Tuple{GeneName, SeqName, Char}

typealias BufOut BufferedStreams.BufferedOutputStream

immutable GeneInfo
   gene::GeneName
   name::SeqName
   strand::Bool # pos is true, neg is false
end

GeneInfo{S <: AbstractString}( gene::S, name::S, strand::Char ) = 
                               GeneInfo( convert(GeneName, gene),
                                         convert(SeqName, name), 
                                         strand == '+' ? true : 
                                                         false )

immutable TxInfo
   name::RefSeqId
   txstart::CoordInt
   txend::CoordInt
   exnum::CoordInt
end


#=
Basic IO Functions
=#

fixpath( str::String ) = abspath( expanduser( str ) )

function isgzipped( filename::String )
   restr = "\.gz\$"
   re = Base.match(Regex(restr), filename)
   return re == nothing ? false : true
end
