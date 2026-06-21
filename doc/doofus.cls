%%======================================================
%% doofus.cls — Custom LaTeX class
%%======================================================

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{doofus}[Custom class for code documentation (DOcumentation OFten USed)]

%% Base class
\LoadClass[a4paper,onecolumn,11pt]{report}

%%======================================================
%% Page geometry
%%======================================================
\RequirePackage[
  left=2cm,
  right=2cm,
  top=2.5cm,
  bottom=3cm,
  includehead,
  includefoot
]{geometry}


%%======================================================
%% Package setup
%%======================================================
\RequirePackage{amsmath}
\RequirePackage{amssymb}
\RequirePackage{listings}
\RequirePackage{courier}
\RequirePackage{xcolor}
\RequirePackage{graphicx}
\RequirePackage{subfigure}
\RequirePackage{bm}
\RequirePackage{hyperref}
\RequirePackage{bookmark}

%%======================================================
%% Hyperref setup
%%======================================================
\hypersetup{
  colorlinks=true,
  linkcolor=black,
  citecolor=black,
  urlcolor=black,
  pdfborderstyle={/S/U/W 0.5},
  pdfborder={0 0 1},
}

%%======================================================
%% Colors — restricted palette for readability
%%======================================================
% Code background (for listings only)
\definecolor{backcolour}{rgb}{0.95,0.95,0.92}
% Grey shades
\definecolor{doofusgray}{gray}{0.3}
\definecolor{doofusgraylight}{gray}{0.6}
% Blue shades
\definecolor{doofusblue}{RGB}{30,90,160}
\definecolor{doofusbluelight}{RGB}{100,160,240}
% Red shades
\definecolor{doofusred}{RGB}{170,40,40}
\definecolor{doofusredlight}{RGB}{220,100,100}

%%======================================================
%% Item symbols
%%======================================================
\renewcommand{\labelitemi}{$\triangleright$}  % Level 1
\renewcommand{\labelitemii}{$\circ$}          % Level 2
\renewcommand{\labelitemiii}{$\bullet$}       % Level 3
\renewcommand{\labelitemiv}{$\ast$}           % Level 4

%%======================================================
%% Default listings configuration
%%======================================================
\lstset{
  basicstyle=\small\ttfamily,
  columns=fullflexible,
  breaklines=true,
  frame=none,
  framesep=3pt,
  framerule=0pt,
}

%%======================================================
%% Custom listing style: mystyle
%%======================================================
\lstdefinestyle{mystyle}{
  basicstyle=\ttfamily\footnotesize,
  columns=fullflexible,
  breaklines=true,
  frame=none,
  framesep=3pt,
  framerule=0pt,
  commentstyle=\color{doofusgray}\itshape,
  keywordstyle=\color{doofusred}\bfseries,
  numberstyle=\tiny\color{doofusgraylight},
  stringstyle=\color{doofusredlight},
  breakatwhitespace=false,
  captionpos=b,
  keepspaces=true,
  numbers=left,
  numbersep=5pt,
  showspaces=false,
  showstringspaces=false,
  showtabs=false,
  tabsize=2
}

%%------------------------------------------------------
%% C++ syntax highlighting style for listings
%% Uses restricted color palette: greys, blues, reds
%%------------------------------------------------------
\lstdefinestyle{cppstyle}{
  language=C++,
  basicstyle=\ttfamily\footnotesize,
  keywordstyle=\color{doofusred}\bfseries,
  commentstyle=\color{doofusgray}\itshape,
  stringstyle=\color{doofusredlight},
  numberstyle=\tiny\color{doofusgraylight},
  numbers=left,
  stepnumber=1,
  numbersep=5pt,
  showstringspaces=false,
  breaklines=true,
  frame=none,
  framerule=0pt,
  framesep=3pt,
  tabsize=2,
  columns=fullflexible
}

%%------------------------------------------------------
%% Title page customization
%%------------------------------------------------------
\newcommand{\doofustitle}[1]{\def\@doofustitle{#1}}
\newcommand{\doofussubtitle}[1]{\def\@doofussubtitle{#1}}
\newcommand{\doofusauthor}[1]{\def\@doofusauthor{#1}}
\newcommand{\doofusdate}[1]{\def\@doofusdate{#1}}
\newcommand{\doofuslogo}[1]{\def\@doofuslogo{#1}}

%% Default values
\doofustitle{Untitled Report}
\doofussubtitle{}
\doofusauthor{Anonymous}
\doofusdate{\today}
\doofuslogo{}

%%------------------------------------------------------
%% Title page definition
%%------------------------------------------------------
\newcommand{\makedoofustitle}{%
  \begin{titlepage}
    \centering
    \vspace*{3cm}

    {\Huge\bfseries \@doofustitle \par}
    \vspace{0.8cm}
    {\Large \@doofussubtitle \par}

    \vspace{2cm}
    {\Large \@doofusauthor \par}

    \vspace{1cm}
    {\large \@doofusdate \par}

    \vspace{2cm}
    \ifx\@doofuslogo\@empty
    \else
      \includegraphics[width=0.3\textwidth]{\@doofuslogo}
    \fi

    \vfill
  \end{titlepage}
  
  % Reset counters so first section starts at 1 (no chapters)
    \setcounter{chapter}{0}
    \renewcommand{\thesection}{\arabic{section}}
    \setcounter{section}{0}
}

%%------------------------------------------------------
%% Syntax presentation commands
%%------------------------------------------------------
%% \syntax{...}: monospace text with semantic coloring, no background.
%% Uses \parbox to allow line wrapping without the colored box.
\newlength{\@syntaxw}
\newcommand{\syntax}[1]{%
  \setlength{\@syntaxw}{\dimexpr\linewidth-2pt\relax}%
  \leavevmode%
  \parbox{\@syntaxw}{\ttfamily\small\raggedright #1}%
  \endgraf\noindent%
}

% Syntax component commands
\newcommand{\return}[1]{{\color{doofusblue}#1}}
\newcommand{\keyword}[1]{{\color{doofusred}\textbf{#1}}}
\newcommand{\type}[1]{{\color{doofusbluelight}#1}}
\newcommand{\param}[1]{{\color{doofusgray}\textit{\textbf{#1}}}}

%%======================================================
%% End of class
%%======================================================
\endinput
