library(RSQLite)

db =  '~/.slurmscripts_debug.sqlite3';
message('Starting')
while( !  file.exists(db ) ) {
	Sys.sleep(10)
}
message('reading database')
con <- dbConnect(RSQLite::SQLite(), db)


## check if sbatch is working:
## first create a siomple script:
OK = FALSE
try({
write( file="simple.sh", c('#!/bin/bash', 'sleep 1'))
OK = system('sbatch simple.sh', intern=T)
OK = length( grep('Submitted batch job', OK) == 1 )
})
OK

runCMD <- function(x) {
	id = x[1]
	cmd = x[2]
	if ( length(grep('^ *sbatch', cmd )) == 1 ) {
		if ( OK ) {
			pid= system( cmd, intern=T)
			pid = as.numeric(stringr::str_replace( pid, 'Submitted batch job ', ''))
			#browser()
			str = paste(sep='', 'UPDATE SlurmScripts SET PID =\'',pid,'\' where id = \'', id,'\'' ) 
			print(str)
			dbSendQuery(con, str )
			return(pid)
		}
		else {
			cmd = stringr::str_replace( cmd, 'sbatch', 'bash')
			system( cmd, intern=T)
			str = paste(sep='', 'UPDATE SlurmScripts SET PID =\'-1\', Finished = 1 where id = \'', id,'\'' ) 
			dbSendQuery(con, str )
		}
	}
	else {
		str = paste(sep='', 'UPDATE SlurmScripts SET PID =\'-100\', Finished = 1 where id = \'', id,'\'' ) 
		dbSendQuery(con, str )
	}
	FALSE
}

updatePIDs <- function( pids) {
	if ( OK ){

	}
	tab = system('squeue', intern=T)
	running = as.numeric(stringr::str_replace_all(unlist(stringr::str_extract_all( tab , '^ *\\d+')), ' +',''))
	for ( finished in setdiff( pids, running )){
		str = paste(sep='', 'UPDATE SlurmScripts SET Finished = 1 where PID = \'', finished,'\'' ) 
		print(str)
		dbSendQuery(con, str )
	}
	running
}
		
pids = NULL

message('Server starting')

while( 1 ) {
	try({
	res <- dbSendQuery(con, 'SELECT id, cmd FROM SlurmScripts WHERE  Finished = 0')
	cmds = dbFetch(res)
	pids =  c ( pids , apply( cmds, 1, runCMD ) )
	dbClearResult(res)
	pids = pids[which(pids!=FALSE)]
	#browser()
	if ( length( pids) > 0 ) {
		## check if they are still active - otherwise update the database
		pids = updatePIDs ( pids )
	}
	Sys.sleep(2)
	dbClearResult(res)
	})
}
if ( file.exists("simple.sh" )){
	unlink("simple.sh")
}
