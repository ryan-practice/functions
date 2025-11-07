search_restriction_sites <- function(sequence, enzymes, use_iupac = TRUE) {
  if (!is.character(sequence) || length(sequence) != 1L) {
    stop("`sequence` must be a single character string.")
  }
  if (is.null(names(enzymes)) || any(names(enzymes) == "")) {
    stop("`enzymes` must be a named character vector: names = enzyme names, values = motifs.")
  }

  seq <- DNAString(toupper(gsub("\\s+", "", sequence)))
  fixed_flag <- !use_iupac

  results <- lapply(names(enzymes), function(en) {
    motif <- toupper(enzymes[[en]])
    pat   <- DNAString(motif)

    # Forward strand matches
    m_f <- matchPattern(pat, seq, fixed = fixed_flag)
    df_f <- if (length(m_f)) {
      data.frame(
        enzyme       = en,
        recog_motif  = motif,
        strand       = "+",
        start        = start(m_f),
        end          = end(m_f),
        match_seq    = as.character(Views(seq, m_f)),
        stringsAsFactors = FALSE
      )
    } else NULL

    # Reverse strand matches (search RC motif on forward sequence)
    pat_rc <- reverseComplement(pat)
    m_r <- matchPattern(pat_rc, seq, fixed = fixed_flag)
    df_r <- if (length(m_r)) {
      data.frame(
        enzyme       = en,
        recog_motif  = motif,   # report motif in forward orientation
        strand       = "-",
        start        = start(m_r),
        end          = end(m_r),
        match_seq    = as.character(Views(seq, m_r)),  # this is the RC match on the forward string
        stringsAsFactors = FALSE
      )
    } else NULL

    rbind(df_f, df_r)
  })

  out <- do.call(rbind, results)
  if (is.null(out)) {
    out <- data.frame(enzyme=character(), recog_motif=character(), strand=character(),
                      start=integer(), end=integer(), match_seq=character(),
                      stringsAsFactors = FALSE)
  } else {
    out <- out[order(out$start, out$end, out$enzyme, out$strand), ]
    rownames(out) <- NULL
  }
  out
}
