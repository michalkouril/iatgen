requireNamespace("stringr")

#' Data analysis function: Processes and cleans raw IAT data
#' @description Prior to running, please see \code{combineIATfourblocks()}. This function processes, cleans, and scores the combined IAT data. In addition, it returns diagnostics (see examples, below). By default, the function implements the D-score algorithm (Greenwald et al., 2003, p 214, center column). Because it assumes users were forced to correct errors, no error penalty is imposed (unless the user requests it; see below). The function can be easily configured to do other scoring procedures as well. The function accepts as an input four vectors of IAT responses (see \code{prac1}, \code{crit1}, \code{prac2}, and \code{crit2}, below). It returns a list containing a variety of IAT variables, including matrices of clean latencies and other information (see below). The most important is \code{clean$D}, which is the final D scores for the analysis. Users can also extract clean block means for each participant using \code{clean$clean.means.prac1}, \code{clean$clean.means.crit1}, \code{clean$clean.means.prac2}, and \code{clean$clean.means.crit2}. Users can extract matrices of clean latencies using \code{clean$clean.latencies.prac1}, \code{clean$clean.latencies.crit1}, etc. Raw latencies can be requested with \code{clean$raw.latencies.prac1}, etc. Users can request to know whether a trial was correct with \code{clean$clean.correct.prac1}, etc. and precisely which stimulus was used on a given trial with \code{clean$clean.stim.number.prac1}, etc. (Stimuli are numbered based on their order entered within each category and following the sequence "positive, negative, tgtA, tgtB". For example, stimulus 1 is the first positive stimulus). See below for more information on what is returned from this function. The data cleaning function adheres to Greenwald et al. (2003; see also Lane et al., 2005, p. 92 for a simplified table of data cleaning steps). There are four main data cleaning options. First, long responses are usually dealt with by setting \code{timeout.drop=TRUE} (enabled by default), which drops individual trials over a given threshold (\code{timeout.ms}, which is 10000 ms by default). Next, overly short responses (i.e., button mashing) are dealt with by setting \code{fastprt.drop=TRUE} (enabled by default), which drops participants who have too many fast responses (more than a \code{fastprt.percent} proportion [default = .10] of responses faster than \code{fastprt.ms} [default = 300 ms]). Alternatively, one can remove individual fast trials by setting \code{fasttrial.drop=TRUE} (disabled by default), which uses a default threshold of \code{fasttrial.ms=400} ms. (This is seldom used but enables users to use alternative scoring methods [e.g., Greenwald et al., 2003, p 214, right column]). Finally, an error penalty is imposed on incorrect responses in some variants. If the IAT forces participants to correct errors, then no error penalty should be imposed (\code{error.penalty=FALSE}, the default setting). However, if participants are not forced to correct errors, one is added. Most common is a 600 ms penalty above the clean block mean (Greenwald et al., 2003), which is done by setting \code{error.penalty.ms=600}, sometimes known as the D600 scoring procedure. Greenwald et al. (2003) also suggested one could use two standard deviations instead of 600 ms, which is done by setting \code{error.penalty.ms="2SD"}. Finally, the function ensures that the data are not corrupted (i.e., JavaScript malfunction on participant's computer when completing the survey) by requiring that only appropriate characters (numbers, commas, "C", "X", and "END) are in the raw data.
#' @param prac1 A vector of one kind of practice responses (e.g., compatible practice), one per participant.
#' @param crit1 A vector of that same kind of critical responses (e.g., compatible critical), one per participant.
#' @param prac2 A vector of the other kind of practice responses (e.g., incompatible practice), one per participant.
#' @param crit2 A vector of that same kind of critical responses (e.g., incompatible critical), one per participant.
#' @param timeout.drop (Required, set \code{TRUE} by default). Tells the procedure to drop trials over a certain duration; recommended by Greenwald et al. (2003).
#' @param timeout.ms (Required if \code{timeout.drop=TRUE}; set to 10000 by default). Following the Greenwald et al. (2003), individual trials over 10000 ms are dropped (scored as missing). Ignored if \code{timeout.drop=FALSE}.
#' @param fastprt.drop (Required, set \code{TRUE} by default). If enabled, follows Greenwald et al. (2003) in which participants who have more than 10 percent of responses (\code{fastprt.percent = .10}) faster than 300 ms (\code{fastprt.ms=300}) are dropped entirely.
#' @param fastprt.percent (Required if \code{fastprt.drop=TRUE}; set to .10 by default). Set the proportion threshold for \code{fastprt.drop}, above. Ignored if \code{fastprt.drop=FALSE}.
#' @param fastprt.ms (Required if \code{fastprt.drop=TRUE}; set to 300 ms by default). Sets the time threshold for for \code{fastprt.drop}, above. Ignored if \code{fastprt.drop=FALSE}.
#' @param fasttrial.drop (Required, set \code{FALSE} by default). Tells the procedure to drop trials under a certain duration. Not recommended but was validated by Greenwald et al. (2003) as an alternative to dropping fast participants.
#' @param fasttrial.ms (Required if \code{fasttrial.drop=TRUE}; set to 400 ms by default). The threshold for \code{fastprt.drop}, above. Ignored if \code{fastprt.drop=FALSE}.
#' @param error.penalty (Required, set \code{FALSE} by default). Logical value stating whether an error penalty is added. This should be disabled if forced error correction was used in the IAT and enabled otherwise (Greenwald et al., 2003).
#' @param error.penalty.ms (Required if \code{error.penalty=TRUE}; set to \code{error.penalty.ms=600} by default). Following the D600 procedure, IAT errors are scored as the correct-trial block mean plus an error penalty of 600 ms. Can be manually set to any desired value. One can also use the 2SD penalty [Greenwald et al., 2003, p 214, right column] by setting  \code{error.penalty.ms="2SD"}. Ignored if \code{error.penalty=FALSE}.
#' @return Returns a list containing several important elements.
#' \code{skipped} is a vector indicating whether the participant completed the IAT or skipped it. They are dropped from analysis if the IAT was skipped.
#' \code{raw.latencies.prac1} is a matrix of the raw latencies in the first practice block prior to any data cleaning.
#' \code{raw.latencies.crit1} is a matrix of the raw latencies in the first critical block prior to any data cleaning.
#' \code{raw.latencies.prac2} is a matrix of the raw latencies in the second practice block prior to any data cleaning.
#' \code{raw.latencies.crit2} is a matrix of the raw latencies in the second critical block prior to any data cleaning.
#' \code{raw.stim.number.prac1} is a matrix of the raw stimuli ID numbers in the first practice block prior to any data cleaning.
#' \code{raw.stim.number.crit1} is a matrix of the raw stimuli ID numbers in the first critical block prior to any data cleaning.
#' \code{raw.stim.number.prac2} is a matrix of the raw stimuli ID numbers in the second practice block prior to any data cleaning.
#' \code{raw.stim.number.crit2} is a matrix of the raw stimuli ID numbers in the second critical block prior to any data cleaning.
#' \code{raw.correct.prac1} is a matrix stating whether each trial was correct (logical) in the first practice block prior to cleaning.
#' \code{raw.correct.crit1} is a matrix stating whether each trial was correct (logical) in the first critical block prior to cleaning.
#' \code{raw.correct.prac2} is a matrix stating whether each trial was correct (logical) in the second practice block prior to cleaning.
#' \code{raw.correct.crit2} is a matrix stating whether each trial was correct (logical) in the second critical block prior to cleaning.
#' \code{timeout.drop} is the logical value stating whether this feature was enabled in the function call (see above).
#' \code{timeout.ms} is the timeout threshold specified in the function call (see above), used if timeout.drop is enabled.
#' \code{num.timeout.removed} is the grand total number of trials removed because they exceeded the timeout threshold in timeout.ms.
#' \code{timeout.rate} is a vector indicating the proportion of responses per participant that were scored as missing due to timeouts.
#' \code{num.timeout.removed.prac1} is the number of trials removed in the first practice block because they exceeded the timeout threshold in timeout.ms.
#' \code{num.timeout.removed.crit1} is the number of trials removed in the first critical block because they exceeded the timeout threshold in timeout.ms.
#' \code{num.timeout.removed.prac2} is the number of trials removed in the second practice block because they exceeded the timeout threshold in timeout.ms.
#' \code{num.timeout.removed.crit2} is the number of trials removed in the second critical block because they exceeded the timeout threshold in timeout.ms.
#' \code{fasttrial.drop} is the logical value stating whether this feature was enabled in the function call (see above).
#' \code{fasttrial.ms} is the time threshold specified in the function call (see above), used if fasttrial.drop is enabled.
#' \code{num.fasttrial.removed} is the grand total number of trials removed because they exceeded the fasttrial threshold in fasttrial.ms.
#' \code{fasttrial.rate} is a vector indicating the percentage of responses per participant that were scored as missing due to rapid speeds.
#' \code{num.fasttrial.removed.prac1} is the number of trials removed in the first practice block because they exceeded the fasttrial threshold in fasttrial.ms.
#' \code{num.fasttrial.removed.crit1} is the number of trials removed in the first critical block because they exceeded the fasttrial threshold in fasttrial.ms.
#' \code{num.fasttrial.removed.prac2} is the number of trials removed in the second practice block because they exceeded the fasttrial threshold in fasttrial.ms.
#' \code{num.fasttrial.removed.crit2} is the number of trials removed in the second critical block because they exceeded the fasttrial threshold in fasttrial.ms.
#' \code{fastprt.drop} is the logical value as specified by the function call (see above).
#' \code{fastprt.ms} is the threshold for as  specified in the function call (see above), used if fastprt.drop is enabled.
#' \code{fastprt.percent} is the proportion of trials specified in the function call (see above).
#' \code{drop.participant} is a logical vector indicating whether the participant's responses have been dropped due to excessive fast responses (if fastprt.drop is enabled).
#' \code{fastprt.count} is the number of participants dropped for excessive fast responding (if fastprt.drop is enabled).
#' \code{fastprt.rate} is the proportion of participants dropped for excessive fast responding (if fastprt.drop is enabled).
#' \code{error.penalty} is a logical value stating whether an error penalty is enabled.
#' \code{error.num.prt} is a vector of the number of erroneous trials per participant (after data cleaning is complete).
#' \code{error.rate.prt} is a vector of the proportion of erroneous trials per participant (after data cleaning is complete).
#' \code{error.rate} is the proportion of the entire set of clean trials which are erroneous trials.
#' \code{error.rate.prac1} is the proportion of the prac1 block set of clean trials which are erroneous trials.
#' \code{error.rate.crit1} is the proportion of the crit1 block set of clean trials which are erroneous trials.
#' \code{error.rate.prac2} is the proportion of the prac2 block set of clean trials which are erroneous trials.
#' \code{error.rate.crit2} is the proportion of the crit2 block set of clean trials which are erroneous trials.
#' \code{clean.latencies.prac1} is a matrix of the clean latencies in the first practice block.
#' \code{clean.latencies.crit1} is a matrix of the clean latencies in the first critical block.
#' \code{clean.latencies.prac2} is a matrix of the clean latencies in the second practice block.
#' \code{clean.latencies.crit2} is a matrix of the clean latencies in the second critical block.
#' \code{clean.stim.number.prac1} is a matrix of the clean stimuli ID numbers in the first practice block.
#' \code{clean.stim.number.crit1} is a matrix of the clean stimuli ID numbers in the first critical block.
#' \code{clean.stim.number.prac2} is a matrix of the clean stimuli ID numbers in the second practice block.
#' \code{clean.stim.number.crit2} is a matrix of the clean stimuli ID numbers in the second critical block.
#' \code{clean.correct.prac1} is a matrix stating whether each trial was correct (logical) in the first practice block.
#' \code{clean.correct.crit1} is a matrix stating whether each trial was correct (logical) in the first critical block.
#' \code{clean.correct.prac2} is a matrix stating whether each trial was correct (logical) in the second practice block.
#' \code{clean.correct.crit2} is a matrix stating whether each trial was correct (logical) in the second critical block.
#' \code{clean.means.prac1} is a vector of clean block mean of latencies in the first practice block, one per participant.
#' \code{clean.means.crit1} is a vector of clean block mean of latencies in the first critical block, one per participant.
#' \code{clean.means.prac2} is a vector of clean block mean of latencies in the second practice block, one per participant.
#' \code{clean.means.crit2} is a vector of clean block mean of latencies in the second critical block, one per participant.
#' \code{diff.prac} is a vector (one per person) of the difference between mean latencies compatible and incompatible (practice) blocks.
#' \code{diff.crit} is a vector (one per person) of the difference between mean latencies compatible and incompatible (critical) blocks.
#' \code{inclusive.sd.prac} is a vector (one per person) of the inclusive SD for the practice trials, per Greenwald et al. (2003).
#' \code{inclusive.sd.crit} is a vector (one per person) of the inclusive SD for the critical trials, per Greenwald et al. (2003).
#' \code{D} is a vector (one per person) of the final D scores (i.e., IAT scores).
#' @references Greenwald, A. G., McGhee, D. E., & Schwartz, J. L. K. (1998). Measuring individual differences in implicit cognition: The Implicit Association Test. \emph{Journal of Personality and Social Psychology, 74}, 1464–1480. https://doi.org/10.1037/0022-3514.74.6.1464
#' @references Greenwald, A. G., Nosek, B. A., & Banaji, M. R. (2003). Understanding and using the Implicit Association Test: I. An improved scoring algorithm. \emph{Journal of Personality and Social Psychology, 85}, 197–216. https://doi.org/10.1037/0022-3514.85.2.197
#' @references Lane, K. A., Banaji, M. R., Nosek, B. A., & Greenwald, A. G. (2007). Understanding and using the Implicit Association Test: IV: What we know (so far) about the method. In B. Wittenbrink & N. Schwarz (Eds.), \emph{Implicit measures of attitudes}. (pp. 59–102). New York, NY: Guilford Press.
#' @references Nosek, B. A., Greenwald, A. G., & Banaji, M. R. (2005). Understanding and using the implicit association test: II. Method variables and construct validity. \emph{Personality and Social Psychology Bulletin, 31}, 166–180. https://doi.org/10.1177/0146167204271418
#' @examples \dontrun{
#'
#' ### CLEAN THE IAT USING THE BUILT IN ERROR PENALTY FOR FORCED-ERROR CORRECTION ###
#' clean <- cleanIAT(dat$compatible.prac, dat$compatible.crit, dat$incompatible.prac, dat$incompatible.crit)
#'
#' ### CLEAN THE IAT USING THE D600 PROCEDURE ###
#' clean <- cleanIAT(dat$compatible.prac, dat$compatible.crit, dat$incompatible.prac, dat$incompatible.crit, error.penalty=TRUE, error.penalty.ms=600)
#'
#' ### CLEAN THE IAT USING THE D2SD PROCEDURE###
#' clean <- cleanIAT(dat$compatible.prac, dat$compatible.crit, dat$incompatible.prac, dat$incompatible.crit, error.penalty=TRUE, error.penalty.ms = "2SD")
#'
#' ### CLEAN THE IAT USING THE D2SD PROCEDURE WITH TRIALS UNDER 400 MS DROPPED ###
#' clean <- cleanIAT(dat$compatible.prac, dat$compatible.crit, dat$incompatible.prac, dat$incompatible.crit, fastprt.drop=FALSE, fasttrial.drop=TRUE, fasttrial.ms=400, error.penalty=TRUE, error.penalty.ms = "2SD")
#'
#' ### EXAMINE CLEAN IAT SCORES
#' clean$D
#'
#' ### EXAMINE IAT DIAGNOSTICS ###
#' # TIMEOUT DROP RATE (% of TRIALS) #
#' clean$timeout.rate
#'
#' # LOWER TAIL DROP RATE (% of TRIALS) - NOTE: DISABLED BY DEFAULT #
#' clean$fasttrial.rate
#'
#' # FAST PARTICIPANT DROP COUNT AND RATE (% of SAMPLE) #
#' clean$fastprt.count
#' clean$fastprt.rate
#'
#' # ERROR RATE #
#' clean$error.rate
#' }


cleanIATnew <- function(prac1, crit1, prac2, crit2, timeout.drop=TRUE, timeout.ms=10000, fasttrial.drop=FALSE, fasttrial.ms=400, fastprt.drop=TRUE, fastprt.percent=.10, fastprt.ms=300, error.penalty=FALSE, error.penalty.ms=600, inclusive.sd=TRUE) {

  if (is.null(prac1)){stop("One of your input variables does not exist. Please check your data / variable names and try again.")}
  if (is.null(prac2)){stop("One of your input variables does not exist. Please check your data / variable names and try again.")}
  if (is.null(crit1)){stop("One of your input variables does not exist. Please check your data / variable names and try again.")}
  if (is.null(crit2)){stop("One of your input variables does not exist. Please check your data / variable names and try again.")}

  if (all(is.na(prac1))){stop("One of your input variables is empty")}
  if (all(is.na(prac2))){stop("One of your input variables is empty")}
  if (all(is.na(crit1))){stop("One of your input variables is empty")}
  if (all(is.na(crit2))){stop("One of your input variables is empty")}

  if (length(unique(c(length(prac1), length(prac2), length(crit1), length(crit2))))!=1) { stop("all blocks equal number of participants") }

  # remove participants that (a) didn't finish a block, (b) have corrupt data
  skipped_bad_data <- gsub("([0-9]+[CX][0-9]+,)+END,#,","",paste(prac1,"#",prac2,"#",crit1,"#",crit2,"#,",sep=','))!=""
  if (all(skipped_bad_data)) { stop("No well formatted data -- please check your input") }
  prac1[skipped_bad_data] <- ""
  prac2[skipped_bad_data] <- ""
  crit1[skipped_bad_data] <- ""
  crit2[skipped_bad_data] <- ""

  # remove ,END
  prac1 <- gsub(",END$", "", prac1)
  prac2 <- gsub(",END$", "", prac2)
  crit1 <- gsub(",END$", "", crit1)
  crit2 <- gsub(",END$", "", crit2)

  # create a block list with all 4 blocks
  blocks.list <- list(prac1=prac1,crit1=crit1,prac2=prac2,crit2=crit2)

  # split on comma, remove trials over median (of non-skipped), pad with NA
  blocks.list.split <- lapply(blocks.list, function(x) { strsplit(x,',') })
  num.raw.trials.blocks <- lapply(blocks.list.split, function(x) { lapply(x, function(y) { length(y[!is.na(y)]) }) })
  blocks.list.split <- lapply(blocks.list.split,
                              function(x) {
                                m<-median(unlist(lapply(x,function(y){ replace(length(y), length(y)==0, NA) })),na.rm = TRUE)
                                # add NAs for short lines
                                lapply(x,function(y){ y[1:m] })
                              })

  # record stim number, correct indicator and latencies
  # before removing timeouts, fast trials and fast participants
  raw.stim.number.blocks <- lapply(blocks.list.split, function(block) {  lapply(block, function(y){ as.numeric(gsub("[CX][0-9]*","", y))}) })
  raw.correct.blocks <- lapply(blocks.list.split, function(block) {  lapply(block, function(y){ gsub("[0-9]*","", y)}) })
  raw.latencies.blocks <- lapply(blocks.list.split, function(block) {  lapply(block, function(y){ as.numeric(gsub("[0-9]*[CX]","", y))}) })
  num.raw.trials  <- rowSums(sapply(as.data.frame(do.call(cbind, num.raw.trials.blocks)), as.numeric))

  # BEGIN identify fast participants
  fastprt.ms.run <- fastprt.ms
  if (fastprt.drop == FALSE) fastprt.ms.run <- -Inf

  # extract trials <fastprt.ms, count them and remove participants who's % exceeds fastprt.percent
  blocks.list.fastprt <- lapply(blocks.list.split, function(block) { lapply(block, function(y){ y[as.numeric(gsub("[0-9]*[CX]","", y))<fastprt.ms.run] }) })
  num.fastprt.removed.blocks <- lapply(blocks.list.fastprt, function(block) { lapply(block, function(y) { length(na.omit(y)) }) })
  number.fastprt <- rowSums(sapply(as.data.frame(do.call(cbind, num.fastprt.removed.blocks)), as.numeric))
  skipped_fast_prts <- (number.fastprt > (num.raw.trials * fastprt.percent))
  # END identify fast participants

  # BEGIN drop long (timeout) trials (and remove them from furher analysis)
  timeout.ms.run <- timeout.ms
  if (timeout.drop == FALSE) timeout.ms.run <- Inf

  blocks.list.timeout <- lapply(blocks.list.split, function(block) { lapply(block, function(y){ ifelse(as.numeric(gsub("[0-9]*[CX]","", y))>timeout.ms.run,y,NA) }) })
  blocks.list.split <- lapply(blocks.list.split, function(block) { lapply(block, function(y){ ifelse(as.numeric(gsub("[0-9]*[CX]","", y))<=timeout.ms.run,y,NA) }) })

  num.timeout.removed.blocks <- lapply(blocks.list.timeout, function(block) { lapply(block, function(y) { sum( !is.na( y ) ) }) })
  num.timeout.removed <- sum(unlist(num.timeout.removed.blocks))
  # END drop long (timeout) trials (and remove them from furher analysis)

  # BEGIN drop fast trials (and remove them from furher analysis)
  fasttrial.ms.run <- fasttrial.ms
  if (fasttrial.drop == FALSE) fasttrial.ms.run <- -Inf

  blocks.list.fast <- lapply(blocks.list.split, function(block) { lapply(block, function(y){ ifelse(as.numeric(gsub("[0-9]*[CX]","", y))<fasttrial.ms.run,y,NA) }) })
  blocks.list.split <- lapply(blocks.list.split, function(block) { lapply(block, function(y){ ifelse(as.numeric(gsub("[0-9]*[CX]","", y))>=fasttrial.ms.run,y,NA) }) })

  num.fasttrial.removed.blocks <- lapply(blocks.list.fast, function(block) { z<-lapply(block, function(y) { length(na.omit(y)) }); sum(unlist(z)) })
  num.fasttrial.removed <- sum(unlist(num.fasttrial.removed.blocks))
  # END drop fast trials (and remove them from furher analysis)

  # calculate rates of dropping
  timeout.rate <- num.timeout.removed / sum(num.raw.trials)
  fasttrial.rate <- num.fasttrial.removed / sum(num.raw.trials)
  fastprt.count <- sum(skipped_fast_prts, na.rm=T)
  fastprt.rate <- sum(skipped_fast_prts, na.rm=T) / sum(!skipped_bad_data, na.rm=T)

  # all dropped participants are those with bad data and fast ones
  drop.participant <- skipped_bad_data | skipped_fast_prts

  # set all values for the dropped participants to NA
  blocks.list.split <- lapply(blocks.list.split, function(block) { mapply(function(x,y){ if(y) { rep(NA,length(x)) } else { x } },x=block, y=drop.participant,SIMPLIFY=FALSE) })

  # extract clean stim number, correct indicator and latencies
  # data from dropped participants will be set to NA because we removed their trials just above
  clean.stim.number.blocks <- lapply(blocks.list.split, function(block) {  lapply(block, function(y){ as.numeric(gsub("[CX][0-9]*","", y))}) })
  clean.correct.blocks <- lapply(blocks.list.split, function(block) {  lapply(block, function(y){ gsub("[0-9]*","", y)}) })
  clean.latencies.blocks <- lapply(blocks.list.split, function(block) {  lapply(block, function(y){ as.numeric(gsub("[0-9]*[CX]","", y))}) })
  num.clean.trials.blocks <- lapply(clean.latencies.blocks, function(block) { lapply(block, function(y) { length(na.omit(y)) })})
  num.clean.trials <- rowSums(sapply(as.data.frame(do.call(cbind, num.clean.trials.blocks)), as.numeric))

  # skip participants that have no data left after eliminating fast and timeout trials
  skipped_no_data_left <- num.clean.trials==0

  # all dropped participants are those with bad data and fast ones
  drop.participant <- skipped_bad_data | skipped_fast_prts | skipped_no_data_left

  # extract number of errors per block (note: length returns 0 for dropped participants too - need to set NA for dropped participants)
  num.errors.blocks <- lapply(clean.correct.blocks, function(block) { lapply(block, function(y) { length(na.omit(y[y=='X'])) })})
  # set all values for the dropped participants to NA
  num.errors.blocks <- lapply(num.errors.blocks, function(block) { mapply(function(x,y){ if(y) { rep(NA,length(x)) } else { x } },x=block, y=drop.participant, SIMPLIFY=FALSE) })

  error.rate.per_blocks <- mapply(function(block.num.errors, block.num.clean.trials) { sum(na.omit(unlist(block.num.errors)))/sum(na.omit(unlist(block.num.clean.trials))) }, block.num.errors=num.errors.blocks, block.num.clean.trials=num.clean.trials.blocks, SIMPLIFY=FALSE)

  # number of errors and rate per participant
  error.num.prt <- rowSums(sapply(as.data.frame(do.call(cbind, num.errors.blocks)), as.numeric))
  error.rate.prt <- error.num.prt / num.clean.trials

  # overall error rate
  error.rate <- sum(na.omit(error.num.prt)) / sum(na.omit(num.clean.trials))

  # extract correct trials and calculate mean
  clean.correct.latencies.blocks <- lapply(blocks.list.split, function(block) { lapply(block, function(y){ as.numeric(gsub("[0-9]*[CX]","",y[gsub("[0-9]*","", y)=='C'])) }) })
  clean.correct.means.blocks <- lapply(clean.correct.latencies.blocks, function(block) { z<-unlist(lapply(block, function(y) { mean(y,na.rm=TRUE)})); z[is.nan(z)]<-NA; names(z)<-seq_along(z); z })

  # calculate stddev for all participant-blocks and correct trials only
  std.nopenalty <- lapply(clean.latencies.blocks, function(block) { lapply(block, function(y) { sd(y,na.rm=TRUE)}) })
  std.correct.nopenalty <- lapply(clean.correct.latencies.blocks, function(block) { lapply(block, function(y) { sd(y,na.rm=TRUE)}) })

  # add error penalty if enabled and error.penalty.ms set to number of "2SD"
  if(error.penalty==TRUE && is.numeric(error.penalty.ms)){
    # only to the X (errors) replace with block mean+penalty.ms
    clean.latencies.blocks<-mapply(function(block,mean,correct) {
      mapply(function(b,m,c) {
        ifelse(c=='C',b,m+error.penalty.ms)
      },b=block,m=mean,c=correct, SIMPLIFY = FALSE)
    },block=clean.latencies.blocks, mean=clean.correct.means.blocks, correct=clean.correct.blocks, SIMPLIFY = FALSE)
  } else if (error.penalty==TRUE && error.penalty.ms=="2SD"){
    # only to the X (errors) replace with block mean+2*SD
    clean.latencies.blocks<-mapply(function(block,mean,correct,block.std) {
      mapply(function(b,m,c,s) {
        ifelse(c=='C',b,m+2*s)
      },b=block,m=mean,c=correct,s=block.std, SIMPLIFY = FALSE)
    },block=clean.latencies.blocks, mean=clean.correct.means.blocks, correct=clean.correct.blocks, block.std=std.correct.nopenalty, SIMPLIFY = FALSE)
  }

  # calcalate mean latency per block per participant, NaN -> NA
  clean.means.blocks <- lapply(clean.latencies.blocks, function(block) { z<-unlist(lapply(block, function(y) { mean(y,na.rm=TRUE)})); z[is.nan(z)]<-NA; names(z)<-seq_along(z); z })

  # prac sd
  inclusive.trials <- mapply(function(block1,block2) { c(block1,block2)}, block1=clean.latencies.blocks$prac1, block2=clean.latencies.blocks$prac2, SIMPLIFY = FALSE)
  inclusive.sd.prac <- unlist(lapply(inclusive.trials, function(y) { sd(y,na.rm=TRUE) }))

  # crit sd
  inclusive.trials <- mapply(function(block1,block2) { c(block1,block2)}, block1=clean.latencies.blocks$crit1, block2=clean.latencies.blocks$crit2, SIMPLIFY = FALSE)
  inclusive.sd.crit <- unlist(lapply(inclusive.trials, function(y) { sd(y,na.rm=TRUE) }))

  ## Dscore

  diff.prac <- clean.means.blocks$prac2-clean.means.blocks$prac1
  diff.crit <- clean.means.blocks$crit2-clean.means.blocks$crit1

  D.prac <- diff.prac / inclusive.sd.prac
  D.crit <- diff.crit / inclusive.sd.crit
  D <- (D.prac + D.crit) / 2

  # to match previous version
  names(error.num.prt) <- seq_along(error.num.prt)
  names(error.rate.prt) <- seq_along(error.rate.prt)
  names(drop.participant) <- seq_along(drop.participant)
  inclusive.sd.prac[is.na(inclusive.sd.prac)] <- 0/(0-1)
  inclusive.sd.crit[is.na(inclusive.sd.crit)] <- 0/(0-1)

  return(list(
    # block.list.fast=block.list.fast,
    skipped=skipped_bad_data,
    raw.latencies.prac1=as.data.frame(do.call("rbind", raw.latencies.blocks$prac1)),
    raw.latencies.crit1=as.data.frame(do.call("rbind", raw.latencies.blocks$crit1)),
    raw.latencies.prac2=as.data.frame(do.call("rbind", raw.latencies.blocks$prac2)),
    raw.latencies.crit2=as.data.frame(do.call("rbind", raw.latencies.blocks$crit2)),
    raw.stim.number.prac1=as.data.frame(do.call("rbind", raw.stim.number.blocks$prac1)),
    raw.stim.number.crit1=as.data.frame(do.call("rbind", raw.stim.number.blocks$crit1)),
    raw.stim.number.prac2=as.data.frame(do.call("rbind", raw.stim.number.blocks$prac2)),
    raw.stim.number.crit2=as.data.frame(do.call("rbind", raw.stim.number.blocks$crit2)),
    raw.correct.prac1=as.data.frame(do.call("rbind", raw.correct.blocks$prac1)),
    raw.correct.crit1=as.data.frame(do.call("rbind", raw.correct.blocks$crit1)),
    raw.correct.prac2=as.data.frame(do.call("rbind", raw.correct.blocks$prac2)),
    raw.correct.crit2=as.data.frame(do.call("rbind", raw.correct.blocks$crit2)),
    timeout.drop=timeout.drop,
    timeout.ms=timeout.ms,
    num.timeout.removed=num.timeout.removed,
    timeout.rate=timeout.rate,
    num.timeout.removed.prac1=sum(unlist(num.timeout.removed.blocks$prac1)),
    num.timeout.removed.crit1=sum(unlist(num.timeout.removed.blocks$crit1)),
    num.timeout.removed.prac2=sum(unlist(num.timeout.removed.blocks$prac2)),
    num.timeout.removed.crit2=sum(unlist(num.timeout.removed.blocks$crit2)),
    fasttrial.drop=fasttrial.drop,
    fasttrial.ms=fasttrial.ms,
    num.fasttrial.removed=num.fasttrial.removed,
    fasttrial.rate=fasttrial.rate,
    num.fasttrial.removed.prac1=num.fasttrial.removed.blocks$prac1,
    num.fasttrial.removed.crit1=num.fasttrial.removed.blocks$crit1,
    num.fasttrial.removed.prac2=num.fasttrial.removed.blocks$prac2,
    num.fasttrial.removed.crit2=num.fasttrial.removed.blocks$crit2,
    fastprt.drop=fastprt.drop,
    fastprt.ms=fastprt.ms,
    fastprt.percent=fastprt.percent,
    drop.participant=drop.participant,
    fastprt.count=fastprt.count,
    fastprt.rate=fastprt.rate,
    error.penalty=error.penalty,
    error.num.prt=error.num.prt,
    error.rate.prt=error.rate.prt,
    error.rate=error.rate,
    error.rate.prac1=error.rate.per_blocks$prac1,
    error.rate.crit1=error.rate.per_blocks$crit1,
    error.rate.prac2=error.rate.per_blocks$prac2,
    error.rate.crit2=error.rate.per_blocks$crit2,
    clean.latencies.prac1=as.data.frame(do.call("rbind", clean.latencies.blocks$prac1)),
    clean.latencies.crit1=as.data.frame(do.call("rbind", clean.latencies.blocks$crit1)),
    clean.latencies.prac2=as.data.frame(do.call("rbind", clean.latencies.blocks$prac2)),
    clean.latencies.crit2=as.data.frame(do.call("rbind", clean.latencies.blocks$crit2)),
    clean.stim.number.prac1=as.data.frame(do.call("rbind", clean.stim.number.blocks$prac1)),
    clean.stim.number.crit1=as.data.frame(do.call("rbind", clean.stim.number.blocks$crit1)),
    clean.stim.number.prac2=as.data.frame(do.call("rbind", clean.stim.number.blocks$prac2)),
    clean.stim.number.crit2=as.data.frame(do.call("rbind", clean.stim.number.blocks$crit2)),
    clean.correct.prac1=as.data.frame(do.call("rbind", clean.correct.blocks$prac1)),
    clean.correct.crit1=as.data.frame(do.call("rbind", clean.correct.blocks$crit1)),
    clean.correct.prac2=as.data.frame(do.call("rbind", clean.correct.blocks$prac2)),
    clean.correct.crit2=as.data.frame(do.call("rbind", clean.correct.blocks$crit2)),
    clean.means.prac1=clean.means.blocks$prac1,
    clean.means.crit1=clean.means.blocks$crit1,
    clean.means.prac2=clean.means.blocks$prac2,
    clean.means.crit2=clean.means.blocks$crit2,
    diff.prac=diff.prac,
    diff.crit=diff.crit,
    inclulsive.sd.prac=inclusive.sd.prac,
    inclusive.sd.crit=inclusive.sd.crit,
    D=D
  ))
}
