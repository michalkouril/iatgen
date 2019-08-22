iatgen (1.2.6)
==============

iatgen (pronounced “I A T gen”) is an R package and Shiny App that
builds and analyzes Qualtrics surveys that contain IATs (Implicit
Association Tests; Greenwald et al., 1998) following a procedure
developed by Carpenter et al. (2018; preprint available at
<a href="https://psyarxiv.com/hgy3z/" class="uri">https://psyarxiv.com/hgy3z/</a>).

Specifically, Carpenter et al. developed procedures for
“survey-software” IATs. These are IATs constructed out of modified
survey elements that have been edited by adding custom JavaScript and
HTML code. The R “iatgen” package was developed as a tool for
customizing and pasting this code into Qualtrics so that functional IAT
surveys can be rapidly built and analyzed.

#### What Exactly Does iatgen Do?

Several frequently asked questions about survey-software IATs and iatgen
can be found at www.iatgen.wordpress.com. We recommend that people take
a look at that web page if questions are not answered here.

First, iatgen is *not* a software tool for running IATs.

Iatgen does not run IATs. Qualtrics does! Instead, iatgen implements the
survey-software IAT method automatically. That is, it configures code
files (HTML, JavaScript) based on your input, copies that code into a
Qualtrics survey template, and outputs a Qualtrics survey containing
your desired IAT–all in just a few seconds. However, importantly, iatgen
is not a ‘software tool’ for running the IAT. It is simply a method for
implementing the survey-software IAT method. Qualtrics is the software
tool that runs the IAT. It would, in theory, be possible for researchers
to skip iatgen and configure Qualtrics IATs manually…but that would be
laborious and prone to error (e.g., typos when coding).

Another benefit of iatgen is that it provides a suite of data analysis
tools for processing the resulting data.

*Please note that the iatgen R package is licensed only for
non-commercial (e.g., academic) use under a Creative Commons (CC BY-NC
4.0) license. More details are provided under Licsense, below.*

Please note that iatgen is not “official” IAT software; that is, it is
neither produced nor endorsed by the IAT creators. Although we have
painstakingly read the IAT literature and faithfully followed those
procedures with our code, the official IAT software remains any software
endorsed by the IAT’s creators. Although the use of Qualtrics as an IAT
tool has been validated by Carpenter et al. (2018), this procedure and
its code were generated by Carpenter et al. and were not provided or
endorsed by the IAT’s creators.

Getting Started
---------------

The purpose of this tutorial is to walk you through how to use the
`iatgen` package to build and analyze IATs.

#### Installation

iatgen can be installed on your computer using the `devtools` package.
You first need to install this package if you do not have it. In
addition, iatgen will use commands from the `stringr` package, so this
should be installed as well.

    install.packages("devtools")
    install.packages("stringr")

Next, iatgen can be installed using the `install_github()` command from
the `devtools` package:

    devtools::install_github("iatgen/iatgen")

#### Loading iatgen

iatgen can be loaded as normal with `library()`:

    library(iatgen)

#### Getting Help

That the primary functions in iatgen have built-in help documentation.
For example, detailed information on `writeIATfull()` can be obtained
with `?writeIATfull()`.

#### Shiny App

Users who do not wish to use the R package can use our Shiny web app,
which has the same features, at
<a href="https://applibs.shinyapps.io/iatui2/" class="uri">https://applibs.shinyapps.io/iatui2/</a>.

Building an IAT
---------------

A brief vignette and summary of IAT build features is provided. Please
read our methods paper
(<a href="https://psyarxiv.com/hgy3z/" class="uri">https://psyarxiv.com/hgy3z/</a>)
for more information.

#### Words-Only Example

All iatgen IATs run in Qualtrics. To build an IAT in Qualtrics, users
must (1) configure JavaScript and HTML files, (2) copy them into
Qualtrics, and (3) configure the Qualtrics files. This is all done
automatically by iatgen, via the `writeIATfull()` function. An example
looks like this.

    writeIATfull(IATname="flowins",
                 posname="Pleasant", 
                 negname="Unpleasant",
                 Aname="Flowers",
                 Bname="Insects",
                 catType="words",
                 poswords = c("Gentle", "Enjoy", "Heaven", "Cheer", "Happy", "Love", "Friend"),
                 negwords = c("Poison", "Evil", "Gloom", "Damage", "Vomit", "Ugly", "Hurt"),
                 tgtType="words",
                 Awords = c("Orchid", "Tulip", "Rose", "Daffodil", "Daisy", "Lilac", "Lily"),
                 Bwords = c("Wasp", "Flea", "Roach", "Centipede", "Moth", "Bedbug", "Gnat"),
                 
                 #advanced options with recommended IAT settings
                 n=c(20, 20, 20, 40, 40, 20, 40),
                 qsf=TRUE, 
                 note=TRUE,
                 correct.error=TRUE,
                 pause=250, 
                 errorpause=300, #not used if correct.error=TRUE
                 tgtCol="black",
                 catCol="green",
                 norepeat=FALSE
    )

I walk through each argument here. More detailed information can be
obtained with `?writeIATfull()`.

-   First, the user specifies `IATname`, which becomes part of the file
    names for the HTML and JavaScript files that are created. Because
    this is used for file names and folder names, it should be short and
    avoid any special characters.

-   Next, we specify the labels that appear on the screen (in the upper
    corners of the IAT). The positive category is named with the
    `posname` argument; the negative category is named with the
    `negname` argument; Target A is named with `Aname`; and Target B is
    named with `Bname`. Note that, for our purposes, we assume a
    “compatible” association is
    `Target A + Positive, Target B + Negative`. In other words, we
    assume Target A will be more associated with the positive category
    and Target B will be more associated with the negative category. At
    the end of the day a positive IAT *D*-score means that one was
    faster in the `Target A + Positive, Target B + Negative`
    (compatible) configuration.

-   Next, we need to tell R whether either (or both) categories and
    targets use text or image stimuli. This is done with the `catType`
    and `tgtType` arguments, which should be set to either `images` or
    `words`. If we set them to images (see below), we will specify the
    stimuli differently than if we set them to words.

-   Next, we set the stimuli. If words, we use the arguments `poswords`,
    `negwords`, `Awords`, and `Bwords` as appropriate. Those should be
    vectors of words or strings of text. If we use images, we instead
    specify `posimgs`/`negimgs` and/or `Aimgs`/`Bimgs` (as appropriate),
    which would be vectors of image URLs. This is illustrated below, so
    don’t worry about this now. This is all you need to set the basic
    information for the IAT.

-   Next, we enter “advanced” options. These are set by default, but
    it’s good to know what your script is doing, so let’s discuss them.
    First, we set `n=c(20, 20, 20, 40, 40, 20, 40)`. This is the number
    of trials per block, a numeric vector of length seven. This means we
    have 20 trials in the first block, 20 in the second block, 20 in the
    third block, 40 in the fourth block (critical combined block), 40 in
    the fifth block (washout block with direction reversed; see Nosek et
    al., 2005 for rationale for setting this to 40), 20 in the sixth
    block, and 40 in the seventh block.

-   The `qsf=TRUE` argument is set (as it is by default) to tell R to
    build a Qualtrics Survey File `*.QSF`. If you don’t want this, set
    it to `FALSE` and R will instead create the JavaScript and HTML
    files manually for you in the user’s working directory (saved as
    .txt files).

-   The `note=TRUE` includes a note telling users what the keys are
    during the task. Some researchers use this to remind users that the
    task uses the “E” and “I” keys. Some people prefer these not on the
    screen.

-   The `correct.error=TRUE` tells R to write the code such that the
    timer continues until participants enter the correct response. Under
    this common IAT variant, errors will result in an error message on
    the screen that persists until the correct response is entered. If
    set to `FALSE`, then R will record whichever answer the user enters
    and display an error message that flashes on the screen only for as
    long as specified by `errorpause` (by default, `errorpause=300`
    milliseconds).

-   The `pause=250` argument tells R to make the duration between trials
    to be 250 milliseconds, or a quarter of a second.

-   The colors of the text stimuli and labels can be set. Typically, in
    an IAT, they are different for targets and categories (to reduce
    confusion). By default, they are set with `tgtCol="black"` and
    `catCol="green"` but can be set to any CSS color name.

-   The `norepeat=FALSE` option uses a random order of presentation of
    trials within each block. Please note that stimuli are *selected*
    for inclusion in the IAT by randomly sampling without replacement
    from stimuli pools (meaning that stimuli will not be selected more
    than once into a set of trials until ALL stimuli from that category
    have been sampled). However, in terms of the order in which those
    stimuli are displayed, setting this to `TRUE` will keep stimuli in
    the order sampled, meaning that a participant will also not *see* a
    duplicate until all other stimuli from that category have veen
    displayed.

#### Image-Based IATs

Targets, categories, or both can use images. Images should be sized 250
x 250 pixels in PNG format and hosted via the user’s Qualtrics account
(tutorial at
<a href="https://osf.io/ntd97/" class="uri">https://osf.io/ntd97/</a>).

Then, `tgtType` and/or `catType` arguments are set to “images” (as
appropriate), and `poswords`/`negwords` are replaced with
`posimgs`/`negimgs` and/or `Awords`/`Bwords` are replaced with
`Aimgs`/`Bimgs` (as appropriate). The only difference between the word
stimuli vectors and the image vectors is that the image vectors are
vectors of image URLs. For stability reasons and on the basis of our own
testing, we *strongly* recommend users of images only host images on
their *own* Qualtrics accounts and follow the guidelines found in the
tutorial referenced above.

Because URLs are long, we recommend specifying vectors of images URLs in
advance and referencing them in the function call:

    goodjpg <- c("www.website.com/gentle.jpg",
                 "www.website.com/enjoy.jpg",
                 "www.website.com/Heaven.jpg",
                 "www.website.com/Cheer.jpg")

    badjpg <- c("www.website.com/Poison.jpg",
                "www.website.com/Evil.jpg.",
                "www.website.com/Vomit.jpg",
                "www.website.com/Ugly.jpg")

    Ajpg <- c("www.website.com/Orchid.jpg",
                 "www.website.com/Tulip.jpg",
                 "www.website.com/Rose.jpg",
                 "www.website.com/Daisy.jpg")

    Bjpg <- c("www.website.com/Wasp.jpg",
                "www.website.com/Flea.jpg",
                "www.website.com/Moth.jpg",
                "www.website.com/Bedbug.jpg")

    writeIATfull(IATname="flowins",
                 posname="Pleasant", 
                 negname="Unpleasant",
                 Aname="Flowers",
                 Bname="Insects",
                 catType="images",
                 posimgs = goodjpg,
                 negimgs = badjpg,
                 tgtType="images",
                 Aimgs = Ajpg,
                 Bimgs = Bjpg,
                 
                 #advanced options with recommended IAT settings
                 n=c(20, 20, 20, 40, 40, 20, 40),
                 qsf=TRUE, 
                 note=TRUE,
                 correct.error=TRUE,
                 pause=250, 
                 errorpause=300, #not used if correct.error=TRUE
                 tgtCol="black",
                 catCol="green"
    )

#### The Qualtrics Survey

Detailed information about this Qualtrics survey is beyond the scope of
this document and is discussed in depth in the Carpenter et al. (2018)
preprint found at
<a href="https://psyarxiv.com/hgy3z/" class="uri">https://psyarxiv.com/hgy3z/</a>.

Of note, however, is that (1) each IAT block is one question and (2)
there are four permutations of the IAT exist, counterbalancing the
left/right starting position for both Target A and the positive
category. Because each IAT consists of 7 blocks, these occupy 28 survey
questions (7 blocks x 4 permutations). These questions are named using
both the question number (Q1-Q28) and a 3-digit code identifying which
IAT permutation it comes from, based on the starting position of Target
A (RP = Target A starts right, initially paired with positive; RN =
starts left with negative; LP = starts left with positive; LN = starts
left with negative). Thus, “Q9 RN2” is the second block in the IAT where
Target A starts on the right side, initially paired with negative (i.e.,
incompatible block comes first). Researchers should carefully consult
our manuscript prior to use.

Analysis
--------

Once data are collected, iatgen can process the resultant data. Several
data-analysis scripts and a user tutorial are provided via
<a href="https://osf.io/ntd97/" class="uri">https://osf.io/ntd97/</a>.
However, a brief analysis vignette is provided here.

In this vignette, users were *not* asked to correct errors and therefore
the “D600” algorithm is used. Note that you need to know how the IAT was
conducted with respect to errors in order to select the correct analysis
procedure. (This is one great reasons to build using an R script, where
you can save a copy or even share your build script on a repository with
materials, data, code, etc.).

Note too that data must be in the *legacy* Qualtrics CSV format. For
data loaded to R, the row containing detailed question information is
also deleted (per usual, when working with Qualtrics data).

First, the data are loaded:

    #### LOAD THE IATGEN PACKAGE ####
    library(iatgen)

    #### READ YOUR DATA HERE AND SAVE IN R AS "DAT" ####
    dat <- read.csv("IAT Flowers Insects.csv", header=T)

To analyse the IAT, the data must be collapsed into four variables
representing practice/critical versions of the compatible and
incompatible blocks. At present, these are scattered across four
hard-coded permutations of the IAT representing left/right
counterbalancing of the starting positions (naming of these variables is
discussed above and in our manuscript at
<a href="https://psyarxiv.com/hgy3z/" class="uri">https://psyarxiv.com/hgy3z/</a>).
The next step in the analysis is to collapse this down using
`combineIATfourblocks()`:

    ### Collapse  IAT data down ####
    dat$compatible.crit <- combineIATfourblocks(dat$Q4.RP4, dat$Q18.LP4, dat$Q14.RN7, dat$Q28.LN7)
    dat$incompatible.crit <- combineIATfourblocks(dat$Q7.RP7, dat$Q21.LP7, dat$Q11.RN4, dat$Q25.LN4)

    ### Collapse  IAT practice blocks ####
    dat$compatible.prac<- combineIATfourblocks(dat$Q3.RP3, dat$Q17.LP3, dat$Q13.RN6, dat$Q27.LN6)
    dat$incompatible.prac <- combineIATfourblocks(dat$Q6.RP6, dat$Q20.LP6, dat$Q10.RN3, dat$Q24.LN3)

Following this, the researcher runs `cleanIAT()`. In this case, the
researcher is careful to set an `error.penalty=TRUE` and
`error.penalty.ms=600` milliseconds, given that participants were not
forced to correct errors (making this the D600 algorithm; had
participants been forced to correct errors, this would have been
`error.penalty=FALSE`, making it the D-built.in.error.penalty
algorithm). This command is done and the result saved to an object for
further use, typically named `clean`.

    ### Clean the IAT ### 
    clean <- cleanIAT(prac1=dat$compatible.prac, 
                      crit1=dat$compatible.crit, 
                      prac2=dat$incompatible.prac, 
                      crit2=dat$incompatible.crit, 
                      
                      timeout.drop=TRUE, 
                      timeout.ms=10000, 
                      
                      fasttrial.drop=FALSE, 
                      
                      fastprt.drop=TRUE, 
                      fastprt.percent=.10, 
                      fastprt.ms=300, 
                      
                      error.penalty=TRUE, 
                      error.penalty.ms=600)

There are a few things to note in the `cleanIAT()`.

-   First, the first four arguments (`prac1`, `crit1`, `prac2`, and
    `crit2`) represent the practice and critical versions of the
    compatible and incompatible blocks, respectively (see above).

-   In addition, following Greenwald et al. (2003)’s *D*-score
    algorithm, we have set a timeout such that trials above 10,000 *ms*
    are removed (`timeout.drop=TRUE`, `timeout.ms=10000`).

-   We have not set individual fast trials to be removed in the same
    manner (`fasttrial.drop=FALSE`) and instead follow the Greenwald et
    al. (2003) procedure of removing all data from participants who have
    more than 10% of responses under 300 ms (`fastprt.drop=TRUE`,
    `fastprt.percent=.10`, `fastprt.ms=300`). Note, however, that you
    could opt to remove fast *trials* instead of fast participants. In
    our experience, fast participants tend to have very high error rates
    (think: a participant pressing buttons randomly and quickly to skip
    the IAT). Thus, dropping fast participants seems to make sense in
    most contexts.

-   Finally, as noted above, we are adding an error penalty in analysis
    of 600 *ms* because participants were not required to correct errors
    (`error.penalty=TRUE`, `error.penalty.ms=600`). Note that the
    `error.penalty` argument is special in that it’s smart. Greenwald et
    al. (2003) suggested researchers could use two standard deviations
    instead of a fixed penalty. We have added that as an option with
    `error.penalty="2SD"`. If you want to disable it (e.g., if
    participants corrected errors), then set it to
    `error.penalty=FALSE`).

The `clean` object is a list containing *many* things For detailed
information see the built-in help file (`?cleanIAT()`). We focus on a
few here.

First, the number of participants who completed the IAT using
`$skipped`, a logical vector indicating whether each person completed an
IAT or not:

    ### NUMBER OF PARTICIPANTS WHO COMPLETED THE IAT ###
    sum(!clean$skipped)

    ## [1] 201

We see here that 201 people (which was the sample size) submitted a
completed IAT for analysis.

Next, we can see the proportion of trials dropped due to exceeding
10,000 ms (as specified in our function call, above) with
`$timeout.rate`:

    ### TIMEOUT DROP RATE (% of TRIALS) ###
    clean$timeout.rate

    ## [1] 0.001285347

As we see here, it is 1/10 of 1% of trials…a very small amount.

Next, we can get diagnostics on the number of participants dropped due
to overly fast responses with `$fastprt.count` and `$fastprt.rate`:

    ### FAST PARTICIPANT 'BUTTON MASHER' DROP COUNT AND RATE (% of SAMPLE) ###
    clean$fastprt.count

    ## [1] 13

    clean$fastprt.rate

    ## [1] 0.06467662

We see this is 13 participants, or approximately 6% of the sample.

If you wanted to know whether individual participants were dropped or
not, simply request `clean$drop.participant` which returns a logical
vector. This can be used, for example, to inspect those responses in
greater detail.

    clean$drop.participant

    ##     1     2     3     4     5     6     7     8     9    10    11    12 
    ##  TRUE  TRUE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE 
    ##    13    14    15    16    17    18    19    20    21    22    23    24 
    ## FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ##    25    26    27    28    29    30    31    32    33    34    35    36 
    ## FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ##    37    38    39    40    41    42    43    44    45    46    47    48 
    ## FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ##    49    50    51    52    53    54    55    56    57    58    59    60 
    ## FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ##    61    62    63    64    65    66    67    68    69    70    71    72 
    ## FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ##    73    74    75    76    77    78    79    80    81    82    83    84 
    ## FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE 
    ##    85    86    87    88    89    90    91    92    93    94    95    96 
    ## FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ##    97    98    99   100   101   102   103   104   105   106   107   108 
    ## FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ##   109   110   111   112   113   114   115   116   117   118   119   120 
    ##  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE 
    ##   121   122   123   124   125   126   127   128   129   130   131   132 
    ## FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE  TRUE 
    ##   133   134   135   136   137   138   139   140   141   142   143   144 
    ## FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ##   145   146   147   148   149   150   151   152   153   154   155   156 
    ## FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ##   157   158   159   160   161   162   163   164   165   166   167   168 
    ## FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ##   169   170   171   172   173   174   175   176   177   178   179   180 
    ## FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE 
    ##   181   182   183   184   185   186   187   188   189   190   191   192 
    ## FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE 
    ##   193   194   195   196   197   198   199   200   201 
    ## FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE

It’s easy to index in this way. For example, imagine we wanted age
statistics on these people:

    dat$age[clean$drop.participant]

    ##  [1] 27 29 28 19 39 26 23 28 23 23 26 28 28

Next, How accurate were our retained participants? This is given with
`$error.rate`:

    ### ERROR RATE ###
    clean$error.rate

    ## [1] 0.07467388

We are less than 10%, which is considered typical for an IAT. The error
rates of individual participants can also be viewed with:

    clean$error.rate.prt

    ##           1           2           3           4           5           6 
    ##          NA          NA          NA 0.016666667 0.025000000          NA 
    ##           7           8           9          10          11          12 
    ## 0.075000000 0.150000000 0.058333333 0.166666667 0.050000000 0.033333333 
    ##          13          14          15          16          17          18 
    ## 0.288135593 0.033333333 0.025000000 0.133333333 0.091666667 0.041666667 
    ##          19          20          21          22          23          24 
    ## 0.133333333 0.058333333 0.041666667 0.025000000 0.008333333 0.068965517 
    ##          25          26          27          28          29          30 
    ## 0.050000000 0.016666667 0.091666667 0.050000000 0.091666667 0.050000000 
    ##          31          32          33          34          35          36 
    ## 0.058333333 0.025000000 0.066666667 0.141666667 0.116666667 0.016666667 
    ##          37          38          39          40          41          42 
    ## 0.183333333 0.050000000 0.041666667 0.041666667 0.091666667 0.116666667 
    ##          43          44          45          46          47          48 
    ## 0.066666667 0.066666667 0.175000000 0.166666667 0.183333333 0.050000000 
    ##          49          50          51          52          53          54 
    ## 0.041666667 0.066666667 0.168067227 0.066666667 0.050000000 0.041666667 
    ##          55          56          57          58          59          60 
    ## 0.116666667 0.058333333 0.033333333 0.058333333 0.050000000 0.041666667 
    ##          61          62          63          64          65          66 
    ## 0.125000000 0.000000000 0.083333333 0.275000000 0.125000000 0.041666667 
    ##          67          68          69          70          71          72 
    ## 0.075000000 0.058333333 0.033333333 0.058333333 0.100000000 0.000000000 
    ##          73          74          75          76          77          78 
    ## 0.116666667 0.041666667 0.041666667 0.058333333 0.025000000          NA 
    ##          79          80          81          82          83          84 
    ##          NA 0.067226891 0.200000000 0.091666667 0.016666667 0.158333333 
    ##          85          86          87          88          89          90 
    ## 0.158333333 0.008333333 0.100000000 0.116666667 0.041666667 0.216666667 
    ##          91          92          93          94          95          96 
    ## 0.116666667 0.041666667 0.025000000 0.066666667 0.091666667 0.083333333 
    ##          97          98          99         100         101         102 
    ## 0.025000000 0.116666667 0.041666667 0.050000000 0.100000000 0.050420168 
    ##         103         104         105         106         107         108 
    ## 0.033333333 0.133333333 0.058333333 0.008333333 0.008333333 0.033333333 
    ##         109         110         111         112         113         114 
    ##          NA 0.041666667 0.158333333 0.025000000 0.041666667 0.109243697 
    ##         115         116         117         118         119         120 
    ## 0.025210084 0.033333333 0.083333333          NA 0.083333333 0.091666667 
    ##         121         122         123         124         125         126 
    ## 0.058333333 0.042735043 0.075000000 0.066666667 0.141666667 0.075000000 
    ##         127         128         129         130         131         132 
    ## 0.025000000 0.083333333          NA 0.158333333 0.066666667          NA 
    ##         133         134         135         136         137         138 
    ## 0.075000000 0.133333333 0.000000000 0.091666667 0.058333333 0.066666667 
    ##         139         140         141         142         143         144 
    ## 0.066666667 0.058333333 0.208333333 0.075000000 0.141666667 0.008333333 
    ##         145         146         147         148         149         150 
    ## 0.025000000 0.175000000 0.116666667 0.200000000 0.050000000 0.000000000 
    ##         151         152         153         154         155         156 
    ## 0.116666667 0.116666667 0.116666667 0.066666667 0.016666667 0.066666667 
    ##         157         158         159         160         161         162 
    ## 0.085470085 0.041666667 0.016666667 0.041666667 0.041666667 0.025000000 
    ##         163         164         165         166         167         168 
    ## 0.075000000 0.100000000 0.025000000 0.108333333 0.133333333 0.016666667 
    ##         169         170         171         172         173         174 
    ## 0.016666667 0.066666667 0.083333333 0.058333333 0.050000000 0.041666667 
    ##         175         176         177         178         179         180 
    ## 0.058333333 0.058333333          NA 0.100000000          NA 0.025000000 
    ##         181         182         183         184         185         186 
    ## 0.066666667 0.050000000 0.033333333 0.041666667 0.050000000          NA 
    ##         187         188         189         190         191         192 
    ## 0.100000000 0.050000000 0.041666667 0.050000000 0.119658120 0.033333333 
    ##         193         194         195         196         197         198 
    ## 0.100000000 0.194915254 0.025000000 0.066666667 0.066666667 0.116666667 
    ##         199         200         201 
    ## 0.025000000 0.058333333 0.108333333

We can see here that dropped participants have an `NA`. However, if you
wanted to know the error rate of dropped participants, you could re-run
`cleanIAT()` without dropping (e.g., saved as `clean.nodrop`) and then
request the error rates from that
(e.g.,`clean.nodrop$error.rate.prt[clean$drop.participant]`).

You can also view the error rate for each of the four combined blocks
with `clean$error.rate.prac1`, `clean$error.rate.crit1`,
`clean$error.rate.prac2` and `clean$error.rate.crit2`:

    clean$error.rate.prac1

    ## [1] 0.05111821

    clean$error.rate.crit1

    ## [1] 0.05004659

    clean$error.rate.prac2

    ## [1] 0.115016

    clean$error.rate.crit2

    ## [1] 0.09090909

Although not the primary means of analysis for the IAT, you see here
just how many errors there are in the `prac2` and `crit2` blocks–the
incompatible block. There are many other elements in this `clean`
object. Take a look at the help file `?cleanIAT()` to see what you can
get. Or, look at the `names()`:

    names(clean)

    ##  [1] "skipped"                     "raw.latencies.prac1"        
    ##  [3] "raw.latencies.crit1"         "raw.latencies.prac2"        
    ##  [5] "raw.latencies.crit2"         "raw.stim.number.prac1"      
    ##  [7] "raw.stim.number.crit1"       "raw.stim.number.prac2"      
    ##  [9] "raw.stim.number.crit2"       "raw.correct.prac1"          
    ## [11] "raw.correct.crit1"           "raw.correct.prac2"          
    ## [13] "raw.correct.crit2"           "timeout.drop"               
    ## [15] "timeout.ms"                  "num.timeout.removed"        
    ## [17] "timeout.rate"                "num.timeout.removed.prac1"  
    ## [19] "num.timeout.removed.crit1"   "num.timeout.removed.prac2"  
    ## [21] "num.timeout.removed.crit2"   "fasttrial.drop"             
    ## [23] "fasttrial.ms"                "num.fasttrial.removed"      
    ## [25] "fasttrial.rate"              "num.fasttrial.removed.prac1"
    ## [27] "num.fasttrial.removed.crit1" "num.fasttrial.removed.prac2"
    ## [29] "num.fasttrial.removed.crit2" "fastprt.drop"               
    ## [31] "fastprt.ms"                  "fastprt.percent"            
    ## [33] "drop.participant"            "fastprt.count"              
    ## [35] "fastprt.rate"                "error.penalty"              
    ## [37] "error.num.prt"               "error.rate.prt"             
    ## [39] "error.rate"                  "error.rate.prac1"           
    ## [41] "error.rate.crit1"            "error.rate.prac2"           
    ## [43] "error.rate.crit2"            "clean.latencies.prac1"      
    ## [45] "clean.latencies.crit1"       "clean.latencies.prac2"      
    ## [47] "clean.latencies.crit2"       "clean.stim.number.prac1"    
    ## [49] "clean.stim.number.crit1"     "clean.stim.number.prac2"    
    ## [51] "clean.stim.number.crit2"     "clean.correct.prac1"        
    ## [53] "clean.correct.crit1"         "clean.correct.prac2"        
    ## [55] "clean.correct.crit2"         "clean.means.prac1"          
    ## [57] "clean.means.crit1"           "clean.means.prac2"          
    ## [59] "clean.means.crit2"           "diff.prac"                  
    ## [61] "diff.crit"                   "inclulsive.sd.prac"         
    ## [63] "inclusive.sd.crit"           "D"

As you can see, there is a lot here (some of it is simply the arguments
you specified as inputs, which is nice to have in the final object in
case you can’t remember what you specified).

We can estimate internal consistency using a procedure described in the
Carpenter et al. (2016) manuscript and De Houwer and De Bruycker (2007)
using the `IATreliability()` command. This returns a number of things,
including `$reliability` which is a split-half reliability estimate:

    ### RELIABILITY ANALYSIS ###
    IATreliability(clean)$reliability

    ## [1] 0.8058219

We see this IAT is estimated at .80, which is great for an IAT.

The method above is somewhat sophisticated, involving sorting trials
into similar groups, then taking alternating trials, scoring the IAT
separately for each half, correlating them, and using a split-half
correction. Please see the De Houwer and De Bruycker (2007) paper for
more information. One advantage to this method is that it analyzes the
reliability of the IAT *D*-score. In other words, it is actually scoring
the IAT. A variant of Cronbach’s alpha can also be used, which simply
lines up pairs of trials (1st trial, 2nd trial, 3rd trial) from the
incompatible and compatible blocks, takes the difference, and uses those
differencde scores in Conbach’s alpha (see Schnabel, Asendorpf, &
Greenwald, 2008). We have built this into our tool as well.

    IATalpha(clean)$alpha.total

    ## Some items ( V2 V3 V4 V5 V7 V8 V9 V10 V14 V15 V16 V17 V21 V34 V37 ) were negatively correlated with the total scale and 
    ## probably should be reversed.  
    ## To do this, run the function again with the 'check.keys=TRUE' optionSome items ( V1.1 V12.1 V26 V27 V33 V39 ) were negatively correlated with the total scale and 
    ## probably should be reversed.  
    ## To do this, run the function again with the 'check.keys=TRUE' option

    ##  raw_alpha
    ##   0.822489

We see here using alpha that we have a score of .82. This is very
similar to the split-half estimate (see our paper for other examples;
they tend to produce highly similar results).

Next, we can examine the scores. The IAT scores are stored as `$D`. It
is common to put them back into one’s datafile, but they can also be
saved and exported to other software (e.g., SPSS; they will line up with
the rows of the source datafile.) A positive score indicates one had a
preference for the compatible block:

    # place back into dat
    dat$D <- clean$D

    # test for IAT effect
    mean(clean$D, na.rm=T)

    ## [1] 0.6137453

    sd(clean$D, na.rm=T)

    ## [1] 0.3622714

    t.test(clean$D)

    ## 
    ##  One Sample t-test
    ## 
    ## data:  clean$D
    ## t = 23.229, df = 187, p-value < 2.2e-16
    ## alternative hypothesis: true mean is not equal to 0
    ## 95 percent confidence interval:
    ##  0.5616230 0.6658675
    ## sample estimates:
    ## mean of x 
    ## 0.6137453

    #cohen d
    mean(clean$D, na.rm=T) / sd(clean$D, na.rm=T)

    ## [1] 1.694159

Here we see that the mean IAT score was *M* = 0.61, *SD* = 0.36,
*t*(187) = 23.23, *p* &lt; .001, 95% CI \[0.56, 0.67\], *d* = 1.69. This
represents a rather large implicit preference for flowers over insects.

There’s much we can do at this point. For example, we could make a
density plot, which shows us that the distribution is fairly symmetrical
but centered well above zero (indicating that the ‘compatible’ block was
indeed easier for participants).

    library(ggplot2)

    ## Registered S3 methods overwritten by 'ggplot2':
    ##   method         from 
    ##   [.quosures     rlang
    ##   c.quosures     rlang
    ##   print.quosures rlang

    ## 
    ## Attaching package: 'ggplot2'

    ## The following objects are masked from 'package:psych':
    ## 
    ##     %+%, alpha

    ggplot(dat, aes(x=D))+
      geom_density(color="black", fill="light blue")+
      theme_light()

    ## Warning: Removed 13 rows containing non-finite values (stat_density).

![](iat_DEMO_files/figure-markdown_strict/unnamed-chunk-21-1.png)

At this point, these responses could be exported to excel with a command
such as `write.csv()` and pasted into another program such as SPSS
(NOTE: user may need to delete `NA` text in missing responses prior to
pasting into SPSS):

    write.csv(clean$D, "iatOUTPUT.csv")

At this point, these *D*-scores can be correlated with other measures or
otherwise analyzed. If uses wish to report the block means by
participant, these can be found as well:

    ### RT DESCRIPTIVES BY BLOCK
    mean(clean$clean.means.crit1, na.rm=T)

    ## [1] 882.2707

    mean(clean$clean.means.crit2, na.rm=T)

    ## [1] 1065.307

    mean(clean$clean.means.prac1, na.rm=T)

    ## [1] 899.6229

    mean(clean$clean.means.prac2, na.rm=T)

    ## [1] 1165.644

    sd(clean$clean.means.crit1, na.rm=T)

    ## [1] 230.9385

    sd(clean$clean.means.crit2, na.rm=T)

    ## [1] 253.1042

    sd(clean$clean.means.prac1, na.rm=T)

    ## [1] 221.7494

    sd(clean$clean.means.prac2, na.rm=T)

    ## [1] 293.3292

License
-------

The iatgen R package (and associated Shiny App) is licensed only for
*non-commercial (e.g., academic) use* under a
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative
Commons Attribution-NonCommercial 4.0 International License</a> (CC
BY-NC 4.0). The tool was created with the intent of making the IAT
accessible to academic researchers who use Qualtrics for online
research.

No warranty is offered and we assume no liability of any kind for any
consequences that may result from using iatgen. This tool can be
modified and distributed with attribution to us, but cannot be used for
commercial purposes. More details are given in the full text of the
license.

Although we believe the IAT can be validly run via Qualtrics (e.g., as
set up via iatgen) and the use of Qualtrics as an IAT tool has been
validated by Carpenter et al. (2018), this procedure and its code is not
provided or endorsed by the IAT’s creators, and all code for this
project was generated by iatgen’s creators. The official IAT sofware
remains any software endorsed by the IAT’s creators. We hold no
copyright to the IAT itself. We are extremely grateful to the IAT’s
creators, especially Tony Greenwald, for inspiring a cohort of young
scientists such as ourselves to study implicit biases and understand why
people do, think, and feel what they do.

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This
work is licensed under a
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative
Commons Attribution-NonCommercial 4.0 International License</a>.

Authors
-------

The iatgen package was built and maintained by Tom Carpenter
(<a href="mailto:tcarpenter@spu.edu" class="email">tcarpenter@spu.edu</a>),
Michal Kouril, Ruth Pogacar, and Chris Pullig. An early prototype of the
HTML and JavaScript were built by Aleksandr Chakroff. Questions
regarding iatgen should be directed to Tom Carpenter.

Acknowledgments
---------------

We would like to express our profound gratitude to Tony Greenwald and
all other IAT scholars who have come before for inspiring our interest
in this project. We also thank Jordan LaBouff and Stephen Aguilar for
contributing validation data and to Naomi Isenberg for help setting up
our website and user tutorials.
