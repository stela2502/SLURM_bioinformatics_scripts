package stefans_libs::database::slurmscripts;


#  Copyright (C) 2010 Stefan Lang

#  This program is free software; you can redistribute it
#  and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation;
#  either version 3 of the License, or (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/>.

use stefans_libs::database::variable_table;
use base variable_table;
use File::Spec::Functions 'catfile';
use DBD::SQLite;
use File::HomeDir;

##use some_other_table_class;

use strict;
use warnings;


sub new {

    my ( $class, $dbh, $debug ) = @_;
    
	
    Carp::confess ("This class creates it's own dbh!")
	  if ( $dbh );

    my ($self, $dbfile);
	if ( $debug ) {
		$dbfile =catfile( File::HomeDir->my_home, '.slurmscripts_debug.sqlite3' );
	}else{
		$dbfile =catfile( File::HomeDir->my_home, '.slurmscripts.sqlite3' );
	}
	$dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","");

    $self = {
		connection => { driver => 'SQLite'},
        debug => $debug,
        dbh   => $dbh
    };

    bless $self, $class if ( $class eq "stefans_libs::database::slurmscripts" );
    $self->init_tableStructure();

    return $self;

}

sub AddDataset{
	my ( $self, $dataset) = @_;
	my $stm = "INSERT INTO SlurmScripts (cmd, Finished) values ( '".$dataset->{'cmd'}."', ".$dataset->{'Finished'}.');' ;
	return $self->{'dbh'}->prepare($stm)->execute();
}

sub  init_tableStructure {
     my ($self, $dataset) = @_;
     my $hash;
     $hash->{'INDICES'}   = [];
     $hash->{'UNIQUES'}   = [];
     $hash->{'variables'} = [];
     $hash->{'table_name'} = "SlurmScripts";

     $self->{'table_definition'} = $hash;
     $self->{'UNIQUE_KEY'} = [ 'id' ];

     push ( @{$hash->{'variables'}},  {
               'name'         => 'creation_time',
               'type'         => 'DATETIME',
               'NULL'         => '1',
				'default' => 'CURRENT_TIMESTAMP',
               'description'  => '',
               'hidden'       => 1
          }
     );

     push ( @{$hash->{'variables'}},  {
               'name'         => 'cmd',
               'type'         => 'BLOB',
               'NULL'         => '0',
               'description'  => '',
               'hidden'       => 0
          }
     );

     push ( @{$hash->{'variables'}},  {
               'name'         => 'PID',
               'type'         => 'INTEGER',
               'NULL'         => '1',
               'description'  => '',
               'hidden'       => 0
          }
     );

	push ( @{$hash->{'variables'}},  {
               'name'         => 'Finished',
               'type'         => 'CHAR',
               'NULL'         => '1',
               'description'  => '',
				'default'     => 0,
               'hidden'       => 0
          }
     );
	
     $self->{'table_definition'} = $hash;

     $self->{'_tableName'} = $hash->{'table_name'}  if ( defined  $hash->{'table_name'} ); # that is helpful, if you want to use this class without any variable tables

     ##now we need to check if the table already exists. remove that for the variable tables!
     unless ( $self->tableExists( $self->TableName() ) ) {
     	$self->create();
     }
     ## Table classes, that are linked to this class have to be added as 'data_handler',
     ## both in the variable definition and here to the 'data_handler' hash.
     ## take care, that you use the same key for both entries, that the right data_handler can be identified.
     #$self->{'data_handler'}->{''} = some_other_table_class->new( );
     return $dataset;
}



sub waiting {
	my $self = shift;
	return $self->get_data_table_4_search(
		{
			'search_columns' =>
			  [ ref($self) . '.Finished' ],
			'where' => [ [ ref($self) . '.Finished', '=', 'my_value' ] ],
		},
		0
	)->Rows();
}


#sub DO_ADDITIONAL_DATASET_CHECKS {
#	my ( $self, $dataset ) = @_;
#
#	$self->{'error'} .= ref($self) . "::DO_ADDITIONAL_DATASET_CHECKS \n"
#	  unless (1);
#
#	return 0 if ( $self->{'error'} =~ m/\w/ );
#	return 1;
#}
#
#sub INSERT_INTO_DOWNSTREAM_TABLES {
#	my ( $self, $dataset ) = @_;
#	 .= '';
#	return 1;
#}
#
#sub post_INSERT_INTO_DOWNSTREAM_TABLES {
#	my ( $self, $id, $dataset ) = @_;
#	$self->{'error'} .= '';
#	return 1;
#}
#
#sub CHECK_BEFORE_UPDATE{
#	my ( $self, $dataset ) = @_;
#
#	$self->{'error'} .= ref($self) . "::DO_ADDITIONAL_DATASET_CHECKS \n"
#	  unless (1);
#
#	return 0 if ( $self->{'error'} =~ m/\w/ );
#	return 1;
#}


sub expected_dbh_type {
	return 'dbh';
	#return 'database_name';
}


1;
