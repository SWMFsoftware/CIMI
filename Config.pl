#!/usr/bin/perl

# Allow in-place editing                                                        
$^I = "";

# Add local directory to search                                                 
push @INC, ".";

use strict;

our @Arguments       = @ARGV;
our $MakefileDefOrig = "src/Makefile.def";
our $Component = "IM";
our $Code      = "CIMI";
our $Show;
our $Help;
our $ERROR;
our $WARNING;
our $Install;
our $Compiler;

our %Remaining;   # Arguments not handled by share/Scripts/Config.pl

# Planet variables
my $ConfigLog = "config.log";
my $Planet;
my $NewPlanet;
my $Grid;
my $NewGrid;

# Grid size variables in src/ModGrid.f90
my $GridFile = "src/ModGrid.f90";
my $GhostCell;
my $nP;
my $nT;
my $nM;
my $nK;
my $NewGhostCell;
my $NewnP;
my $NewnT;
my $NewnM;
my $NewnK;
my $ShowGhostCell;
my $ShowGrid;

$NewPlanet="EarthHO" if $Install;

# Make sure that config.log exists
`touch $ConfigLog`;

#CIMI non-SWMF developers should use the following gitdir
my $GITDIR   = "https://github.com/MSTEM-QUDA";
#SWMF developers should use the following gitdir
#my $GITDIR   = "git\@github.com:SWMFsoftware";

my $config   = "share/Scripts/Config.pl";
my $gitclone = "share/Scripts/gitclone -s";

# Git clone missing directories as needed. Start with share/ to get $gitclone.
if (not -f $config and not -f "../../$config"){
    `git clone $GITDIR/share; git clone $GITDIR/util`;
}


if(-f $config){
    require $config;
}else{
    require "../../$config";
}

&print_help if $Help;

# Select planet: Current choices are Earth with  H+,O+,He+,e or H+,O+,e, or H+,e
# Select Grid:   Current choices are default and expanded
# More planet configurations will be added in future
foreach (@Arguments){
    if(/^-earthhohe$/i)         {$NewPlanet="EarthHOHe";           next};
    if(/^-earthswho$/i)         {$NewPlanet="EarthSwHO";           next};
    if(/^-earthho$/i)           {$NewPlanet="EarthHO";             next};
    if(/^-earthh$/i)            {$NewPlanet="EarthH";              next};
    if(/^-griddefault$/i)       {$NewGrid="GridDefault";           next};
    if(/^-gridexpanded$/i)      {$NewGrid="GridExpanded";          next};
    if(/^-gridsouthexpanded$/i) {$NewGrid="SouthGridExpanded";     next};
    if(/^-griduniforml$/i)      {$NewGrid="GridUniformL";          next};
    if(/^-ng$/i)                {$ShowGhostCell=1;                 next};
    if(/^-ng=(.*)$/i)           {$NewGhostCell=$1;                 next};
    if(/^-np$/i)                {$ShowGrid=1;                      next};
    if(/^-nt$/i)                {$ShowGrid=1;                      next};
    if(/^-nm$/i)                {$ShowGrid=1;                      next};
    if(/^-nk$/i)                {$ShowGrid=1;                      next};
    if(/^-np=(.*)$/i)           {$NewnP=$1;                        next};
    if(/^-nt=(.*)$/i)           {$NewnT=$1;                        next};
    if(/^-nm=(.*)$/i)           {$NewnM=$1;                        next};
    if(/^-nk=(.*)$/i)           {$NewnK=$1;                        next};
    if(/^-s$/)                  {$Show=1;                          next};
    if(/^-compiler=(.*)$/i)     {$Compiler=$1;                     next};        
    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}

&get_settings;

&set_planet if $NewPlanet and $NewPlanet ne $Planet;

&set_grid   if $NewGrid and $NewGrid ne $Grid;

# Set grid size after set_grid, which overwrites src/ModGrid.f90
&get_grid_size;
&set_grid_size if length($NewGhostCell) or length($NewnP) or length($NewnT)
    or length($NewnM) or length($NewnK);

# &set_compiler  if $Install; # SAMI3 is not present

print "nGhostCell=$GhostCell\n" if $ShowGhostCell;
print "np=$nP,nt=$nT,nm=$nM,nk=$nK\n" if $ShowGrid;

&show_settings if $Show;

exit 0;
############################################################################

sub get_settings{

    open(FILE, $ConfigLog) or 
	die "$ERROR could not open $ConfigLog\n";
    
    while(<FILE>){
        if(/\s*PLANET\s*=\s*(\w+)/){
	    $Planet=$1;
	    last;
	}
    }
    close FILE;
    
    open(FILE, $ConfigLog) or 
	die "$ERROR could not open $ConfigLog\n";
    
    while(<FILE>){
	if(/\s*Grid\s*=\s*(\w+)/){
	    $Grid=$1;
	    last;
	}
    }
    close FILE;

    #die "$ERROR could not find PLANET name in $ConfigLog\n"
    #	unless $Planet;
}

#############################################################################

sub get_grid_size{

    # Read nGhostCell, np, nt, nm and nk from src/ModGrid.f90
    $GhostCell = undef;
    $nP        = undef;
    $nT        = undef;
    $nM        = undef;
    $nK        = undef;

    open(FILE, $GridFile) or die "$ERROR could not open $GridFile\n";
    while(<FILE>){
	next if /^\s*!/;
	s/!.*//;
	$GhostCell = $1 if /\bnGhostCell\s*=\s*(\d+)/i;
	$nP        = $1 if /\bnp\s*=\s*(\d+)/i;
	$nT        = $1 if /\bnt\s*=\s*(\d+)/i;
	$nM        = $1 if /\bnm\s*=\s*(\d+)/i;
	$nK        = $1 if /\bnk\s*=\s*(\d+)/i;
    }
    close FILE;
}

#############################################################################

sub set_grid_size{

    die "$ERROR -ng=$NewGhostCell must be an integer >= 2\n"
	if length($NewGhostCell)
	and ($NewGhostCell !~ /^\d+$/ or $NewGhostCell < 2);
    die "$ERROR -np=$NewnP must be a positive integer\n"
	if length($NewnP) and ($NewnP !~ /^\d+$/ or $NewnP < 1);
    die "$ERROR -nt=$NewnT must be a positive integer\n"
	if length($NewnT) and ($NewnT !~ /^\d+$/ or $NewnT < 1);
    die "$ERROR -nm=$NewnM must be a positive integer\n"
	if length($NewnM) and ($NewnM !~ /^\d+$/ or $NewnM < 1);
    die "$ERROR -nk=$NewnK must be a positive integer\n"
	if length($NewnK) and ($NewnK !~ /^\d+$/ or $NewnK < 1);

    $GhostCell = $NewGhostCell if length($NewGhostCell);
    $nP        = $NewnP        if length($NewnP);
    $nT        = $NewnT        if length($NewnT);
    $nM        = $NewnM        if length($NewnM);
    $nK        = $NewnK        if length($NewnK);

    print "Writing nGhostCell=$GhostCell\n" if length($NewGhostCell);
    print "Writing grid np,nt,nm,nk=$nP,$nT,$nM,$nK\n"
	if length($NewnP) or length($NewnT) or length($NewnM) or length($NewnK);

    my ($FoundGhostCell, $FoundNp, $FoundNt, $FoundNm, $FoundNk);
    @ARGV = ($GridFile);
    while(<>){
	if(not /^\s*!/){
	    # Only modify the code part, not the trailing comment
	    my ($Code, $Comment) = /^([^!]*)(.*)$/s;
	    $FoundGhostCell = 1 if length($NewGhostCell)
		and $Code =~ s/\b(nGhostCell\s*=\s*)\d+/$1$GhostCell/i;
	    $FoundNp = 1 if length($NewnP)
		and $Code =~ s/\b(np\s*=\s*)\d+/$1$nP/i;
	    $FoundNt = 1 if length($NewnT)
		and $Code =~ s/\b(nt\s*=\s*)\d+/$1$nT/i;
	    $FoundNm = 1 if length($NewnM)
		and $Code =~ s/\b(nm\s*=\s*)\d+/$1$nM/i;
	    $FoundNk = 1 if length($NewnK)
		and $Code =~ s/\b(nk\s*=\s*)\d+/$1$nK/i;
	    $_ = $Code.$Comment;
	}
	print;
    }

    die "$ERROR could not find nGhostCell in $GridFile\n"
	if length($NewGhostCell) and not $FoundGhostCell;
    die "$ERROR could not find np in $GridFile\n"
	if length($NewnP) and not $FoundNp;
    die "$ERROR could not find nt in $GridFile\n"
	if length($NewnT) and not $FoundNt;
    die "$ERROR could not find nm in $GridFile\n"
	if length($NewnM) and not $FoundNm;
    die "$ERROR could not find nk in $GridFile\n"
	if length($NewnK) and not $FoundNk;
}

#############################################################################

sub set_planet{
    my $File;
    $Planet = $NewPlanet;

    my $Dir = "src";
    die "$ERROR Directory $Dir is missing\n" unless -d $Dir;

    $File .= " $Dir/ModEarthSwHO.f90" if $Planet eq "EarthSwHO";
    $File .= " $Dir/ModEarthHOHe.f90" if $Planet eq "EarthHOHe";
    $File .= " $Dir/ModEarthHO.f90"   if $Planet eq "EarthHO";
    $File .= " $Dir/ModEarthH.f90"    if $Planet eq "EarthH";

    &shell_command("cp $File src/ModPlanet.f90");

    &shell_command("echo PLANET=$NewPlanet > config.log");
    &shell_command("echo Grid=$NewGrid >> config.log");
}
#############################################################################

sub set_compiler{

    my $Dir = "srcSAMI3";
    die "$ERROR Directory $Dir is missing\n" unless -d $Dir;
    my $File = "$Dir/Makefile.suff";

    if(-e "$File\.$Compiler"){
	&shell_command("cp $File\.$Compiler $File");
    }else{
	&shell_command("rm -f $File; touch $File");
    }
}
#############################################################################


sub set_grid{
    my $File;
    $Grid = $NewGrid;

    my $Dir = "src";
    die "$ERROR Directory $Dir is missing\n" unless -d $Dir;

    $File = " $Dir/ModGrid_default.f90"   if $Grid eq "GridDefault";
    $File = " $Dir/ModGrid_expanded.f90"  if $Grid eq "GridExpanded";
    $File = " $Dir/ModGrid_southexpanded.f90"  if $Grid eq "SouthGridExpanded";
    $File = " $Dir/ModGrid_UniformL.f90"  if $Grid eq "GridUniformL";

    &shell_command("cp $File src/ModGrid.f90");

    &shell_command("echo PLANET=$NewPlanet > config.log");
    &shell_command("echo Grid=$NewGrid >> config.log");
}
############################################################################

sub show_settings{
    print "Planet=$Planet\n";
    print "Grid=$Grid\n";
    print "nGhostCell=$GhostCell\n";
    print "np=$nP,nt=$nT,nm=$nM,nk=$nK\n";
}

############################################################################

sub print_help{
    print "Additional options for CIMI/Config.pl:

-EarthSwHO       Configure CIMI for Earth with Solar Wind H+,H+,O+,e. 

-EarthHOHe       Configure CIMI for Earth with H+,O+,He+,e. 

-EarthHO         Configure CIMI for Earth with H+,O+,e. 

-EarthH          Configure CIMI for Earth with H+,e. 

-GridDefault     Configure CIMI with the default grid

-GridExpanded    Configure CIMI with a grid expanded to higher lat


-GridUniformL    Configure CIMI to use uniform L parameter for drift in lat

-SouthGridExpanded    Configure CIMI with S hemisiphere and a  grid expanded to higher lat

-ng              Print the current number of ghost cells (nGhostCell).

-ng=NG           Set nGhostCell=NG in src/ModGrid.f90.

-np, -nt, -nm, -nk
                 Print the current grid size np,nt,nm,nk.

-np=NP           Set np=NP (latitude grid size) in src/ModGrid.f90.
                 Note: the non-uniform-L grids read np+2 values from
                 xlat_data, so changing np also requires updating xlat_data.

-nt=NT           Set nt=NT (local-time grid size) in src/ModGrid.f90.

-nm=NM           Set nm=NM (magnetic moment grid size) in src/ModGrid.f90.

-nk=NK           Set nk=NK (K grid size) in src/ModGrid.f90.

Note: selecting a grid with -Grid... copies a template over src/ModGrid.f90,
which resets nGhostCell, np, nt, nm and nk unless -ng, -np, -nt, -nm and -nk
are passed again.

These flags are case insensitive.

In the future, other ions could be added, or other planets such as Jupiter or 
Saturn could be included.

EXAMPLE: Reconfigure CIMI for EarthH and show selected planet:

    Config.pl -EarthH -s

EXAMPLE: Set the grid size to np=76, nt=96 and show it:

    Config.pl -np=76 -nt=96 -np
";
    exit 0;
}
