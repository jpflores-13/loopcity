.onLoad <- function(libname, pkgname) {
    # Seqinfo 1.0.0 was split from GenomeInfoDb, creating two separate seqinfo
    # generics. InteractionSet registers its GInteractions method against
    # Seqinfo::seqinfo, but mariner imports seqinfo from GenomeInfoDb (a different
    # generic object). Register a bridge method against GenomeInfoDb::seqinfo so
    # that mariner::pullHicMatrices can resolve chromosome information correctly.
    gid_seqinfo <- GenomeInfoDb::seqinfo
    if (!existsMethod(gid_seqinfo, "GInteractions")) {
        setMethod(gid_seqinfo, "GInteractions", function(x) {
            getFromNamespace("seqinfo", "Seqinfo")(InteractionSet::regions(x))
        })
    }
}
