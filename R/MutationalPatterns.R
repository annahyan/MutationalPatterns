# Global package variables

# Default genome

DEFAULT_GENOME = "BSgenome.Hsapiens.UCSC.hg19"

# Default colours for mutation profile plots
COLORS6 = c(
    "#2EBAED", "#000000", "#DE1C14",
    "#D4D2D2", "#ADCC54", "#F0D0CE")

COLORS7 = c(
    "#2EBAED", "#000000", "#DE1C14",
    "#E98C7B", "#D4D2D2", "#ADCC54",
    "#F0D0CE")

COLORS10 = c(
    "#A7CBDD", "#2776AE", "#B1DE8B",
    "#389C2F", "#F49C9B", "#D41F1E",
    "#FABC73", "#FB8101", "#C3AFD5",
    "#7D599E"
)

COLORS_INDEL_PREDEF = c(
    "#F7BF80", "#ED8212", 
    "#B5D988", "#31A12C", 
    "#E44A39", "#B81C20"
)

COLORS_INDEL = c(
    "#F7BF80", "#ED8212", "#B5D988", "#31A12C",
    "#F8CAB9", "#EA8E77", "#E44A39", "#B81C20",
    "#D0E1F2", "#97C1DE", "#4B97CA", "#1C68AA",
    "#D0CFD4", "#B2AEC5", "#8079AE", "#634298"
)

COLORS_INDEL_COSMIC = c(
  "#F7BF80", "#ED8212", "#B5D988", "#31A12C",
  "#F8CAB9", "#EA8E77", "#E44A39", "#B81C20",
  "#D0E1F2", "#97C1DE", "#4B97CA", "#1C68AA",
  "#D0CFD4", "#B2AEC5", "#8079AE", "#634298"
)

# Predefined substitutions
SUBSTITUTIONS = c('C>A','C>G','C>T','T>A','T>C','T>G')
SUBSTITUTIONS_96 = rep(SUBSTITUTIONS, each=16)
SUBSTITUTIONS_192 = rep(SUBSTITUTIONS, each=32)

C_TRIPLETS = c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT")

T_TRIPLETS = c(
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT")

CONTEXTS_96 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))

# combine substitutions and context in one 
TRIPLETS_96 = paste(substr(CONTEXTS_96,1,1), "[", SUBSTITUTIONS_96, "]", substr(CONTEXTS_96,3,3), sep = "")
TRIPLETS_192_trans = paste(rep(TRIPLETS_96, each=2), c("transcribed", "untranscribed"), sep="-")
TRIPLETS_192_rep = paste(rep(TRIPLETS_96, each=2), c("left", "right"), sep="-")

DBS = c(
  'AC>CA', 'AC>CG', 'AC>CT', 'AC>GA', 'AC>GG', 'AC>GT', 'AC>TA', 'AC>TG', 'AC>TT', 
  'AT>CA', 'AT>CC', 'AT>CG', 'AT>GA', 'AT>GC', 'AT>TA', 
  'CC>AA', 'CC>AG', 'CC>AT', 'CC>GA', 'CC>GG', 'CC>GT', 'CC>TA', 'CC>TG', 'CC>TT',
  'CG>AT', 'CG>GC', 'CG>GT', 'CG>TA', 'CG>TC', 'CG>TT', 
  'CT>AA', 'CT>AC', 'CT>AG', 'CT>GA', 'CT>GC', 'CT>GG', 'CT>TA', 'CT>TC', 'CT>TG',
  'GC>AA', 'GC>AG', 'GC>AT', 'GC>CA', 'GC>CG', 'GC>TA', 
  'TA>AT', 'TA>CG', 'TA>CT', 'TA>GC', 'TA>GG', 'TA>GT',
  'TC>AA', 'TC>AG', 'TC>AT', 'TC>CA', 'TC>CG', 'TC>CT', 'TC>GA', 'TC>GG', 'TC>GT', 
  'TG>AA', 'TG>AC', 'TG>AT', 'TG>CA', 'TG>CC', 'TG>CT', 'TG>GA', 'TG>GC', 'TG>GT',
  'TT>AA', 'TT>AC', 'TT>AG', 'TT>CA', 'TT>CC', 'TT>CG', 'TT>GA', 'TT>GC', 'TT>GG')

DBS_strand = rep(DBS, each = 2)
DBS_trans = paste(rep(DBS, each=2), c("transcribed", "untranscribed"), sep="-")
DBS_rep = paste(rep(DBS, each=2), c("left", "right"), sep="-")

SUBSTITUTIONS_DBS = c("AC>NN", "AT>NN", "CC>NN","CG>NN","CT>NN","GC>NN","TA>NN","TC>NN","TG>NN","TT>NN")
ALT_DBS = do.call(rbind, strsplit(DBS, ">"))[,2]

INDEL = "cosmic"

INDEL_CONTEXT_PREDEF = c(
  paste0('del.rep.len.', 1:5),
  paste0('ins.rep.len.', 1:5),
  paste0('del.mh.bimh.', 1:5),
  paste0('ins.mh.bimh.', 1:5),
  paste0('del.none.len.', 1:5),
  paste0('ins.none.len.', 1:5)
)

INDEL_CLASS_PREDEF = c(
  rep('del.rep', 5), rep('ins.rep', 5),
  rep('del.mh', 5), rep('ins.mh', 5),
  rep('del.none', 5), rep('ins.none', 5)
)

INDEL_CLASS_HEADER_PREDEF = NULL

INDEL_CONTEXT = c(
  paste0('del.1bp.homopol.C.len.', c(1:5,"6+")),
  paste0('del.1bp.homopol.T.len.', c(1:5,"6+")),
  paste0('ins.1bp.homopol.C.len.', c(0:4,"5+")),
  paste0('ins.1bp.homopol.T.len.', c(0:4,"5+")),
  paste0('del.rep.len.2.rep.', c(1:5,"6+")),
  paste0('del.rep.len.3.rep.', c(1:5,"6+")),
  paste0('del.rep.len.4.rep.', c(1:5,"6+")),
  paste0('del.rep.len.5+.rep.', c(1:5,"6+")),
  paste0('ins.rep.len.2.rep.', c(0:4,"5+")),
  paste0('ins.rep.len.3.rep.', c(0:4,"5+")),
  paste0('ins.rep.len.4.rep.', c(0:4,"5+")),
  paste0('ins.rep.len.5+.rep.', c(0:4,"5+")),
  paste0('del.mh.len.2.bimh.1'),
  paste0('del.mh.len.3.bimh.', c(1,2)),
  paste0('del.mh.len.4.bimh.', c(1:3)),
  paste0('del.mh.len.5+.bimh.', c(1:4,"5+"))
)

INDEL_CONTEXT_COSMIC = c(
  paste0('del.1bp.homopol.C.len.', c(1:5,"6+")),
  paste0('del.1bp.homopol.T.len.', c(1:5,"6+")),
  paste0('ins.1bp.homopol.C.len.', c(0:4,"5+")),
  paste0('ins.1bp.homopol.T.len.', c(0:4,"5+")),
  paste0('del.rep.len.2.rep.', c(1:5,"6+")),
  paste0('del.rep.len.3.rep.', c(1:5,"6+")),
  paste0('del.rep.len.4.rep.', c(1:5,"6+")),
  paste0('del.rep.len.5+.rep.', c(1:5,"6+")),
  paste0('ins.rep.len.2.rep.', c(0:4,"5+")),
  paste0('ins.rep.len.3.rep.', c(0:4,"5+")),
  paste0('ins.rep.len.4.rep.', c(0:4,"5+")),
  paste0('ins.rep.len.5+.rep.', c(0:4,"5+")),
  paste0('del.mh.len.2.bimh.1'),
  paste0('del.mh.len.3.bimh.', c(1,2)),
  paste0('del.mh.len.4.bimh.', c(1:3)),
  paste0('del.mh.len.5+.bimh.', c(1:4,"5+"))
)

INDEL_CLASS = c(
  rep("C", 6), rep("T", 6), 
  rep("C", 6), rep("T", 6),
  rep("2", 6), rep("3", 6), rep("4", 6), rep("5+", 6),
  rep("2", 6), rep("3", 6), rep("4", 6), rep("5+", 6),
  rep("2", 1), rep("3", 2), rep("4", 3), rep("5+", 5)
)

INDEL_CLASS_HEADER = c(
  rep("del.1bp", 12), rep("ins.1bp", 12),
  rep("del.rep", 24), rep("ins.rep", 24),
  rep("del.mh", 11)
)

INDEL_CLASS_COSMIC = c(
  rep("C", 6), rep("T", 6), 
  rep("C", 6), rep("T", 6),
  rep("2", 6), rep("3", 6), rep("4", 6), rep("5+", 6),
  rep("2", 6), rep("3", 6), rep("4", 6), rep("5+", 6),
  rep("2", 1), rep("3", 2), rep("4", 3), rep("5+", 5)
)

INDEL_CLASS_HEADER_COSMIC = c(
  rep("del.1bp", 12), rep("ins.1bp", 12),
  rep("del.rep", 24), rep("ins.rep", 24),
  rep("del.mh", 11)
)

INDEL_MATRIX = NULL

# Strand information 
STRAND = rep(c("U","T"), 96)
STRAND_DBS = rep(c("U","T"), 78)
DNA_BASES = c("A", "C", "G", "T")
