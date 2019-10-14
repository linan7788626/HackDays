########################################################################
# Install using conda, I did not make it using this way.
########################################################################
conda config --add channels http://eupsforge.net/conda/dev
conda install lsst-sims

source /Users/<username>/anaconda/bin/eups-setups.sh
setup lsst_sims
conda update lsst-sims

########################################################################
# Install from source
########################################################################
-- Create a new user.
-- unset LSST_HOME EUPS_PATH LSST_DEVEL EUPS_PKGROOT REPOSITORY_PATH
-- mkdir -p /the/LSST/installation/root && cd /the/LSST/installation/root
-- curl -OL https://sw.lsstcorp.org/eupspkg/newinstall.sh
-- bash newinstall.sh
-- source $INSTALL_DIR/loadLSST.bash # for bash users
-- eups distrib install -t v10_1 lsst_apps
-- eups distrib install lsst_sims -t sims
-- setup lsst_sims

# Setup environment
source "/Users/lsststack/lsst/loadLSST.bash"  
setup lsst_sims

########################################################################
# Running examples of CatSim_GalSim interface
########################################################################

# building ssh Tunnel
# Do not forget to send Scott your ssh public key
ssh -L 51433:fatboy-private.phys.washington.edu:1433 simsuser@gateway.astro.washington.edu

# editing config/db.py
vi $SIMS_CATUTILS_DIR/config/db.py
#--------------------------------------------------
root.driver='mssql+pymssql'
root.host='localhost'
root.port='51433'
root.database='LSST'
#--------------------------------------------------

# Edit username and password
vi .lsst/db-auth.paf 
#--------------------------------------------------
database: {  
    authInfo: {
        host: localhost
        port: 51433
        user: LSST-2
        password: L$$TUser
    }
}
#--------------------------------------------------
# Make it safe
cd $HOME
chmod 700 .lsst
chmod 600 .lsst/db-auth.paf

# Refs
https://confluence.lsstcorp.org/display/SIM/Accessing+the+UW+CATSIM+Database
