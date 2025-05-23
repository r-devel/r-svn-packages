% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/cluster/inst/NEWS.Rd"))!
\name{NEWS}
\title{News for \R Package \pkg{cluster}}% MM: look into ../svn-log-from.all
\encoding{UTF-8}

%% NB: The date (yyyy-mm-dd) is the "Packaged:" date in ../DESCRIPTION

%% NB 2: do not release within 1 year before we have bug fix or new feature.

      %% FIXME ?:  pam() is 3 x  SLOWER than in 2021 only on nb-mm5  ?!?!?
      % nb-mm5:  ~/R/Pkgs/cluster.Rcheck-64b.bak/tests/pam.Rout
      % already  cluster 2.1.7 is  10--25 % slower  (but that's "nothing"
      % compared with the slow down to 2.1.9 which uses dys_2(.)

\section{Changes in version 2.1.9 (2025-03-11, svn r8499)}{
  \subsection{Tweaks in C}{
    \itemize{
      \item Remove more \code{<var>--} in C code; try improving speed of C
      level \code{dysta()}.

      \item Using new \code{dys_2(dys, i,j)} macro to avoid accessing
      \sQuote{dys[0] :== 0}, so we got rid of it (see \dQuote{prepend '0.'}
      below).
    }
  }
  \subsection{Misc}{
    \itemize{
      \item use complete argument names in most cases to satisfy lovers of
      \code{options(warnPartialMatch* = TRUE)}.

      \item \code{agnes()}, \code{diana()}, and \code{pam()} and C-level
      \code{twins()}, \code{dysta()}, etc, no longer prepend \code{0.} to
      the distances (which is an allocation waste for largish \code{n}).

      \item updated messages (for translation) from source.
    }
  }
}

\section{Changes in version 2.1.8 (2024-12-10, svn r8468)}{
  \subsection{Bug Fixes}{
     \itemize{
       \item C codes underlying \code{agnes()} and \code{daisy()} under
       \command{UBSAN} and \command{LTO} no longer show pointer expression
       overflow (of a length 0 vector \code{a} with \code{a--} call).
     }
  }
}

\section{Changes in version 2.1.7 (2024-12-06, svn r8466)}{
  \subsection{Misc}{
    \itemize{
      \item \file{man/pam.Rd} \code{set.seed()} for \emph{data}, hence
      reproducible example.

      \item \file{src/daisy.c} and \file{src/dysta.c} produced from their
      Fortran *.f versions by \command{f2c}, my perl script
      \command{f2c-clean} and manual pretty editing.

      C function \code{dysta()} now returns the \code{jhalt} "flag"; it is
      no longer callable from \R, nor is \code{dysta3()} which is 'static'
      now, in \file{src/fanny.c}.

      \bold{Goal}: Use \code{.Call()} instead of \code{.C()} and hence allow long
      vectors for the length \eqn{n(n-1)/2} distance vectors.

      \item \file{fanny-ex.R} tests: use IGNORE_RDIFF_BEGIN .. and some
      explicit cheks to become less platform dependent.
    }
  }
}

\section{Changes in version 2.1.6 (2023-11-30, svn r8311)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{clara(*, metric="gower")} is too buggy, and hence
      \bold{removed} for now, see below (2.1.5).
    }
  }
  \subsection{Misc}{
    \itemize{
      \item tweak \file{tests/fanny-ex.\{R,Rout.save\}} for less "diffs".
    }
  }
}

\section{Changes in version 2.1.5 (2023-11-27, svn r8301)}{
  \subsection{New Features}{
    \itemize{
      \item \code{pam()} and \code{clara()} also accept a
      \code{\link{logical}} matrix as input.
      \item \code{clara()} gets new \code{metric = "gower"} donated by
      Kasper Fischer-Rasmussen.
      --- was removed in 2.1.6 because of valgrind-detected memory leaks
      (and problems/bugs found earlier in R experiments).
    }
  }
  \subsection{Misc}{
    \itemize{
      \item renamed internal function \code{silhouette.default.R} to \code{silhouetteR}.
      \item \file{CITATION} switched from old \code{citEntry()} to \code{bibentry()}.
      \item \code{clara(*, trace = <pos>)} prints slightly better info.
      \item \file{DESCRIPTION} gets \kbd{Enhances: } for packages
      cross-referenced in \file{man/*.Rd}.
    }
  }
}

\section{Changes in version 2.1.4 (2022-08-19, svn r8119)}{
  \subsection{New Features}{
    \itemize{
      \item New \code{medoids()} to compute \code{pam}-consistent medoids,
      given a clustering (and data).
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{help(dissimilarity.object)} % ../man/dissimilarity.object.Rd
      now correctly describes the meaning of a \code{method = "gower"}
      result type \code{"T"}.
      Also, \code{help(daisy)} now extensively documents the different
      options for \code{type = *}, notably \code{"ordratio"} and
      \code{"logratio"}.

      \item \code{pam(x, diss=TRUE)} uses \code{inherits(x, "dissimilarity")} instead
      of semi-deprecated \code{data.class(x) == "dissimilarity")} to check.

      \item Using \code{deparse1()} from \R version 4.0.0 and newer (even
      when \pkg{cluster} is installed into an older versions of R),
      e.g., in plot default titles.
    }
  }
}

\section{Changes in version 2.1.3 (2022-03-28, svn r8062)}{
  \subsection{New Features}{
    \itemize{
      \item \code{as.data.frame(<silhouette>)} now works; suggested by
	Jonathan Marshall (@Massey, NZ).
    }
  }
  \subsection{Translations}{
    \itemize{
      \item Updated German translations (Detlef Steuer)
    }
  }
  \subsection{Tweaks}{
    \itemize{
      \item use explicit cast, quieting \code{-Wconversion}.
    }
  }
}

\section{Changes in version 2.1.2 (2021-04-16, svn r7959)}{% CRAN published 2021-04-17; R 4.1.0
  \subsection{New Features}{
    \itemize{
      \item Updated the Italian translations (by Daniele Medri), the
      French (by Philippe Grosjean) and the German (by Detlef Steuer).
      New Lithuanian translations by Gabriele Stupuriene & Rimantas Zakauskas.
    }
  }
}

\section{Changes in version 2.1.1 (2021-02-11, svn r7925)}{
  \subsection{New Features}{
    \itemize{
      \item Added Italian translations from Daniele Medri, and updated the
      French ones from Philippe Grosjean.

      \item \code{clara()} gains a \code{cluster.only} argument, entirely
      \dQuote{in parallel} to \code{pam()}.

      \item the \code{silhouette()} method for \code{clara} now allows
      \code{full} to be a number in [0,1] instead of just a logical and it gains
      gains new optional argument \code{subset} to specify the indices of
      observations for which silhouette values shold be computed.

      \item Faster \code{pam()} - \sQuote{fasterPAM} from Erich Schubert.
      New optional arguments \code{medoids = "random"}, \code{nstart}, and
      \code{variant}, notably \code{variant = "faster"}.
    }
  }
}

\section{Changes in version 2.1.0 (2019-06-14, svn r7675)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{volume(obj)} is now correct also for ellipsoids of
      dimension \eqn{d > 2}.

      \item \command{--enable-lto} compilation revealed missing
      \code{weights} argument in \file{src/cluster.h}'s declaration of
      \code{cldaisy()}.
    }
  }
  \subsection{Tweaks}{
    \itemize{
      \item Replaced many \dQuote{old style} \verb{`<word>'} quotations.
    }
  }
}

\section{Changes in version 2.0.9 (2019-05-01, svn r7663)}{
  \subsection{Tweaks}{
    \itemize{
      \item \file{src/mona.c} now uses correct \emph{void} \code{F77_NAME}
      and source cleaned.
    }
  }
}

\section{Changes in version 2.0.8 (2019-04-02, svn r7643)}{
  \subsection{New Features}{
    \itemize{
      \item \code{pam()} gets new \code{pamonce} options, for values 3,
      4, or 5, choosing versions of \command{fastpam}, contributed by Erich
      Schubert, Univ. Dortmund.

      \item update \file{tests/*} to work with R >= 3.6.0 (sample.kind)
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item correct thinko in \code{?clusGap}, the help page.
    }
  }
}

\section{Changes in version 2.0.7 (2018-03-29, svn r7509)}{
  \subsection{New Features}{
    \itemize{
      \item \code{clara()} gets new option \code{metric = "jaccard"},
      contributed  by Kamil Kozlowski and Kamil Jadszko.
      %% FIXME:  Also add for  pam() !!

      \item \code{pam()} and \code{clara()} use \code{match.arg(metric)}
      and hence \code{metric} can be abbreviated (and invalid strings
      give an error instead of being interpreted as \code{"euclidean"}).
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item The bug fix of \code{clara(*, correct.d = TRUE)} (from
      version 2.0.4) for the NA-data case now also applies to the
      internal C function \code{selec()}.
    }
  }
}

\section{Changes in version 2.0.6 (2017-03-10, svn r7332)}{
  \subsection{New Features}{
    \itemize{
      \item \code{mona()} now C- instead of Fortran-based (having used
      f2c etc) and now has a \code{trace.lev} option which allows
      progress reporting
      \dQuote{remembers} if the original data had missing values.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{mona(<1-column>)} no longer loops infinitely but signals
      an error.
    }
  }
}

\section{Changes in version 2.0.5 (2016-10-07, svn r7278)}{
  \subsection{New Features}{
    \itemize{
      \item \code{clusGap()} gets a new option \code{scaleH0}, and
      \code{scaleH0 = "original"} is an alternative to the default PCA
      rotation.%% still see ../TODO-MM !

      \item \code{clusGap()} now also stores its \code{call} and uses
      that for \code{print()}ing and (by default in the \code{main} title)
      for \code{plot()}ing \code{"clusGap"} objects.

      \item
      __ MOSTLY NOT IMPLEMENTED yet __ %%% TODO !!!

      \code{diana()} gets new optional argument \code{stop.at.k}.
      When a positive integer, the DIANA algorithm will stop early, as
      much desirable for large \eqn{n}.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{daisy()} gets 3+1 new options \code{warn*} which allow
      to suppress three different kind of warnings, as these are
      undesirable in some cases.  With thanks to Kirill Müller for the
      convincing context.

      \item \code{pam()} now signals an error when there are more than
      65536 observational units (whereas it could segfault previously),
      thanks to a patch from Mikko Korpela, Helsinki.
    }
  }
}

\section{Changes in version 2.0.4 (2016-04-16, svn r7186)}{
  \subsection{New Features}{
    \itemize{
      \item \code{clusGap()} gets a new option \code{d.power = 1}
      allowing to choose the basic weight statistic as it was originally
      proposed, namely \emph{squared} distances by setting \code{d.power = 2}.
      %% ~/R/MM/Pkg-ex/cluster/Gonzalez-on-clusGap.R <--
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item fix small glitch in silhouette's help page.

      \item Finally fixed a bug (in the original Fortran code from
      Rousseeuw!) in clara's distance computation when there are
      \code{NA}s in the data.  As the fix is not backward compatible,
      a warning is produced (for the time being) if there \emph{are}
      \code{NA}s and the user does not explicitly use \code{clara(*, correct.d = TRUE)}.
    }
  }
}

\section{Changes in version 2.0.3 (2015-07-20, svn r6985)}{
  \subsection{New Features}{
    \itemize{
      \item This new \file{NEWS.Rd} file -- going to replace \file{ChangeLog}
      eventually.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item import all we need (but not more) from the "base" pkgs
      (stats, graphics, ...).
    }
  }
}

\section{Changes in version 2.0.2 (2015-06-18, svn r6955)}{
  \subsection{New Features}{
    \itemize{
      \item using new \code{anyNA()} where appropriate.
      \item New Korean translations, thanks to Chel Hee Lee.
      \item \code{plotpart()}: \code{cmdscale()} tweaks.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item valgrind detected missing allocation (\code{nisol["1"]} for k=1).
      \item typo R/daisy.q (R bug %once we require R >= 3.2.0: \PR{16430}
      \Sexpr[results=rd]{tools:::Rd_expr_PR(16430)}).
    }
  }
}

\section{Changes in version 2.0.1 (2015-01-31, svn r6877)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item Fix \code{silhouette( obj )} for \code{obj <- pam(x, k = 1)}.
    }
  }
}

\section{Changes in version 2.0.0 (2015-01-29, svn r6874)}{
  \subsection{New Features}{
    \itemize{
      \item \code{pam()} now using \code{.Call()} instead of
      \code{.C()} is potentially considerably more efficient.
      \item \code{agnes()} has improved \code{trace} behaviour; also,
      some invalid \code{par.method = *} settings now give an early and
      understandable error message.
      \item \code{lower.to.upper.tri.inds()} (etc) now returns \code{integer}.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{.C(..)} and \code{.Fortran(..)}: no longer using
      \code{DUP=FALSE} as that has become deprecated.
    }
  }
}

\section{Changes in version 1.15.3 (2014-09-04, svn r6804)}{
  \subsection{New Features}{
    \itemize{
      \item \code{agnes()} and \code{diana()} finally get, respectively
      work with a \code{trace.lev} option.
      \item \code{plot.(agnes|diana)()} now deals well with long
      \code{call}s, by using multiple title lines.
      \item Message translations now also for C level error messages.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{agnes(*, method="flexible", par.method = c(a1, a2, b, c))},
      i.e., \code{length(alpha) == 4}, finally works \emph{correctly}.
    }
  }
}

\section{Changes in version 1.15.2 (2014-03-31, svn r6724)}{
  \subsection{New Features}{
    \itemize{
      \item Rewrote parts of the R level messages so they are more
      easily translatable, thanks to proposals by Lukasz Daniel.
      \item French translations from Philippe Grosjean.
    }
  }
}

\section{Changes in version 1.15.1 (2014-03-13, svn r6676)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{mona} example not working in \R < 3.0.x.
    }
  }
}

\section{Changes in version 1.15.0 (2014-03-11, svn r6672)}{
  \subsection{New Features}{
    \itemize{
      \item \code{agnes(*, method = "gaverage")} contributed by Pierre
      Roudier.
      \item documentation improvements;
      \item better translatable messages and translation updates.
    }
  }
}


%% ============================== FIXME ===========================
%%        ~~~~~~~~~
%% use ../ChangeLog
%%        ~~~~~~~~~
%% and then
%%
%% use ../svn-log-from.all
%%        ~~~~~~~~~~~~~~~~
%% and ../../cluster_Archive.lst  {~= CRAN  src/contrib/Archive/cluster/ :
%%
\section{Changes in version 1.14.4 (2013-03-26, svn r6469)}{
  \subsection{Tweaks}{
    \itemize{
      \item 2 de translations; update \file{mona.Rout.save}
    }
  }
}

\section{Changes in version 1.14.3 (2012-10-14, svn r6318)}{
  \subsection{New Features}{
    \itemize{
      \item Polnish translations from Lukasz Daniel.
    }
  }
  \subsection{Misc}{
    \itemize{
      \item Tweaks for translations.
      \item Improve \code{clusGap()} examples.
    }
  }
}


\section{Changes in version 1.14.2 (2012-02-06, svn r6069)}{
  \subsection{New Features}{
    \itemize{
      \item New \code{clusGap()} to compute the \dQuote{cluster Gap}
      goodness-of-fit statistic.
      \item Export \code{coefHier()}
      \item \code{as.dendrogram()} now for \code{"twins"}, not just
      \code{"agnes"}.
    }
  }
}

\section{Changes in version 1.14.1 (2011-10-16, svn r5928)}{
  \subsection{New Features}{
    \itemize{
      \item First translations (into German, thanks to Detlef Steuer).
      \item better \code{citation("cluster")}
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{plot.silhouette(..., col = <one per cluster>)} had
      ordering bug.
    }
  }
}

\section{Changes in version 1.14.0 (2011-06-07, svn r5842)}{
  \subsection{New Features}{
    \itemize{
      \item Make \code{clara()} more \code{pam()}-like;
      provide \code{pamLike = *} option.
      \item \code{clusplot.default()} re-factored; use new
      \code{mkCheckX()}; new args \code{add = FALSE}, \code{cex},
      \code{cex.txt}.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{agnes( <n=1> )}: correct error now.
    }
  }
}
%% -> /sfs/w/ftp/CRAN/src/contrib/Archive/cluster/

%%  214765  Feb 21  2011   cluster_1.13.3.tar.gz
%%  213663  Nov 10  2010   cluster_1.13.2.tar.gz
%%  214083  Jun 25  2010   cluster_1.13.1.tar.gz
%%  214677  Apr  2  2010   cluster_1.12.3.tar.gz
%%  214577  Oct  6  2009   cluster_1.12.1.tar.gz
%%  215041  May 13  2009   cluster_1.12.0.tar.gz
%%  211085  Mar 31  2009   cluster_1.11.13.tar.gz
%%  321990  Jan  7  2009   cluster_1.11.12.tar.gz
%%  245055  Jun 16  2008   cluster_1.11.11.tar.gz
%%  243446  Feb 29  2008   cluster_1.11.10.tar.gz
%%  216573  Oct  2  2007   cluster_1.11.9.tar.gz
%%  215257  Sep  4  2007   cluster_1.11.8.tar.gz
%%  216815  Jun  5  2007   cluster_1.11.7.tar.gz
%%  216729  Apr 27  2007   cluster_1.11.6.tar.gz
%%  211615  Mar 31  2007   cluster_1.11.5.tar.gz
%%  211634  Dec 12  2006   cluster_1.11.4.tar.gz
%%  203692  Dec  2  2006   cluster_1.11.3.tar.gz
%%  210927  Sep  7  2006   cluster_1.11.2.tar.gz
%%  210091  Aug 25  2006   cluster_1.11.1.tar.gz
%%  210215  May 18  2006   cluster_1.11.0.tar.gz
%%  195962  Mar 21  2006   cluster_1.10.5.tar.gz
%%  197577  Jan 26  2006   cluster_1.10.4.tar.gz
%%  197853  Jan 26  2006   cluster_1.10.3.tar.gz
%%  190839  Aug 31  2005   cluster_1.10.2.tar.gz
%%  190975  Jul  3  2005   cluster_1.10.1.tar.gz
%%  189042  Jun 13  2005   cluster_1.10.0.tar.gz
%%  179723  Apr  4  2005   cluster_1.9.8.tar.gz
%%  176832  Jan 24  2005   cluster_1.9.7.tar.gz
%%  174742  Aug 24  2004   cluster_1.9.6.tar.gz
%%  174218  Aug  4  2004   cluster_1.9.5.tar.gz
%%  175565  Jun 26  2004   cluster_1.9.4.tar.gz
%%  173097  Jun 18  2004   cluster_1.9.3.tar.gz
%%  173251  Jun 13  2004   cluster_1.9.2.tar.gz
%%  169773  Apr 12  2004   cluster_1.9.1.tar.gz
%%  170071  Mar 14  2004   cluster_1.8.1.tar.gz
%%  165322  Jan 22  2004   cluster_1.8.0.tar.gz
%%  161548  Sep 24  2003   cluster_1.7.6.tar.gz
%%  161359  Sep  3  2003   cluster_1.7.5.tar.gz
%%  161257  Jul 18  2003   cluster_1.7.4.tar.gz
%%  160252  Jun 11  2003   cluster_1.7.3.tar.gz
%%  158265  Jun  4  2003   cluster_1.7.2.tar.gz
%%  157386  May  1  2003   cluster_1.7.1.tar.gz
%%  155161  Mar 26  2003   cluster_1.7.0.tar.gz
%%  154089  Dec 31  2002   cluster_1.6-4.tar.gz
%%  154987  Dec  5  2002   cluster_1.6-3.tar.gz
%%  154261  Oct 23  2002   cluster_1.6-2.tar.gz
%%  147063  Sep 10  2002   cluster_1.6-1.tar.gz
%%  131808  Jul 30  2002   cluster_1.5-2.tar.gz
%%  116292  Jun 19  2002   cluster_1.5-1.tar.gz
%%  113972  Mar 31  2002   cluster_1.4-2.tar.gz
%%  113889  Mar  7  2002   cluster_1.4-1.tar.gz
%%  116698  Jan 24  2002   cluster_1.4-0.tar.gz
%%  105552  Dec 19  2001   cluster_1.3-6.tar.gz
%%  105390  Nov  7  2001   cluster_1.3-5.tar.gz
%%  105275  Aug 24  2001   cluster_1.3-4.tar.gz
%%  103626  Jun  8  2001   cluster_1.3-3.tar.gz
%%   99698  Jan  4  2001   cluster_1.3-2.tar.gz
%%   91608  Feb 18  2000   cluster_1.2-3.tar.gz
%%   91736  Dec 29  1999   cluster_1.2-2.tar.gz
%%   93048  Dec  5  1999   cluster_1.2-1.tar.gz

%% ============================== FIXME ===========================



\section{Version 1.2-1 -- 1.13-3}{
  \subsection{Many historical releases}{
    \itemize{
      \item 60 more CRAN releases of the package \pkg{cluster}
      from Dec 1999 to Feb 2011, see also the \file{ChangeLog} file and
      \command{svn log}.
    }
  }
}

% How can I add vertical space ?
% \preformatted{} is not allowed, nor is \cr


\section{Version 1.2-0 (1999-04-11)}{
  \subsection{First CRAN release of the \pkg{cluster} package, by Kurt Hornik}{
    \itemize{
      \item Martin Maechler had its own version independently.
      \item Both closely modeled after \code{clus} the tarball off JSS.
  }}

  \subsection{R Functions -- Fortran Files}{
    \itemize{
      \item \code{agnes()} -- \file{twins.f} for the \dQuote{twins} \code{agnes} and \code{diana}.
      \item \code{clara()} -- \code{clara.f}
      \item \code{daisy()} -- \file{daisy.f} (and \file{meet.f})
      \item \code{diana()} -- (twins.f)
      \item \code{fanny()} -- \file{fanny.f}
      \item \code{mona()}  -- \file{mona.f}
      \item \code{pam()}   -- \file{pam.f}
    }
  }
  \subsection{Data Sets}{
    \itemize{
      \item agriculture
      \item animals
      \item flower
      \item ruspini
      \item votes.repub
    }
  }

  \subsection{Further Features}{
    \itemize{
      \item all Examples in \file{man/*.Rd} hand edited to become
      executable.
      \item \code{summary()}, \code{print()} (and
    \code{print.summary.**()} methods) for the six basic \R functions above.
    }
  }
}



\section{Version 1.1-2 (1998-06-16)}{
  \subsection{Renamed previous \pkg{clus} to \pkg{cluster}}{
    \itemize{ \item . }
  }
}
\section{Version 1.1-1 (1998-06-15)}{
  \subsection{New Features}{
    \itemize{
      \item started \file{ChangeLog}
    }
  }
}
