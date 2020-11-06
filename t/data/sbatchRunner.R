library(RSQlite)

while( !  file.exists( '~/.slurmscripts.sqlite3' ) ) {
	Sys.sleep(10)
}

con <- dbConnect(RSQLite::SQLite(), '~/.slurmscripts.sqlite3')}


## check if sbatch is working:
## first create a siomple script:
write( file="simple.sh", c('#!/bin/bash', 'sleep 1'))
OK = system('sbatch simple.sh', intern=T)
OK = length( grep('Submitted batch job', OK) == 1 )

OK

runCMD <- function(x) {
	id = x[1]
	cmd = x[2]
	if ( length(grep('^ *sbatch', cmd )) == 1 ) {
		if ( OK ) {
			pid= system( cmd, intern=T)
		}
		else {
			cmd = stringr::str_replace( cmd, 'sbatch', 'bash')
			system( cmd, intern=T)
			str = paste(sep='', 'UPDATE SlurmScripts SET PID =\'-1\', Finished = 1 where id = \'', id,'\'' ) 
			dbSendQuery(con, str )
		}
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

unlink("simple.sh")