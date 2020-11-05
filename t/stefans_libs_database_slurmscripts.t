#! /usr/bin/perl
use strict;
use warnings;
use Test::More tests => 5;
BEGIN { use_ok 'stefans_libs::database::slurmscripts' }

my ( $value, @values, $exp );
my $obj = stefans_libs::database::slurmscripts -> new(undef, 1);
is_deeply ( ref($obj) , 'stefans_libs::database::slurmscripts', 'simple test of function stefans_libs::database::slurmscripts -> new()' );

$obj-> create();
$obj-> AddDataset({'cmd' => "Not to be run", 'Finished' => 0 } );

is_deeply ( $obj->waiting(), 1, "One waiting");

$obj-> AddDataset({'cmd' => "Not to be run", 'Finished' => 0 } );
$obj-> AddDataset({'cmd' => "Not to be run", 'Finished' => 0 } );
$obj-> AddDataset({'cmd' => "Not to be run", 'Finished' => 0 } );

is_deeply ( $obj->waiting(), 4, "4 waiting");

$obj-> create();

is_deeply ( $obj->waiting(), 0, "zero waiting");


#print "\$exp = ".root->print_perl_var_def($value ).";\n";


my $example_R_function = "

library(RSQlite)

while( !  file.exists( '~/.slurmscripts.sqlite3' ) ) {
	Sys.sleep(10)
}

con <- dbConnect(RSQLite::SQLite(), '~/.slurmscripts.sqlite3')}


runCMD <- function(x) {
	id = x[1]
	cmd = x[2]
	if ( length(grep('^ *sbatch', cmd )) == 1 ) {
		pid= system( cmd, intern=T)
		if ( length( grep('^Submitted batch job ', pid)) == 1 ) {
			pid = as.numeric(stringr::str_replace( pid, 'Submitted batch job ', ''))
			#browser()
			str = paste(sep='', 'UPDATE SlurmScripts SET PID =\'',pid,'\' where id = \'', id,'\'' ) 
			print(str)
			dbSendQuery(con, str )
			return(pid)
		}
	}
	FALSE
}

updatePIDs <- function( pids) {
	tab = system('squeue', intern=T)
	running = as.numeric(stringr::str_replace_all(unlist(stringr::str_extract_all( tab , '^ *\\d+')), ' +',''))
	for ( finished in setdiff( tab, running )){
		str = paste(sep='', 'UPDATE SlurmScripts SET Finished = 1 where PID = \'', finished,'\'' ) 
		dbSendQuery(con, str )
	}
	running
}
		
pids = NULL


while( 1 ) {
	try({
	res <- dbSendQuery(con, 'SELECT id, cmd FROM SlurmScripts WHERE  Finished = 0')
	cmds = dbFetch(res)
	pids =  c ( pids , apply( cmds, 1, runCMD ) )
	dbClearResult(res)
	pids = pids[which(pids!=FALSE)]
	if ( length( pids) > 0 ) {
		## check if they are still active - otherwise update the database
		pids = updatePIDs ( pids )
	}
	Sys.sleep(2)
	dbClearResult(res)
	})
}

";