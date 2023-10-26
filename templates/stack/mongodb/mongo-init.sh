set -e

UUID=$(cat /proc/sys/kernel/random/uuid)

# Secrets must be read from secret database env file as they
# can be exposed in Docker commands and configurations when
# configured directly in the container deployment

source $MONGODB_CONFIG_FILE

mongosh <<EOF
disableTelemetry()

use admin
db.createUser({
  user:  '$MONGODB_USERNAME',
  pwd: '$MONGODB_PASSWORD',
  roles: [
    {role: 'readWriteAnyDatabase', db: 'admin'},
    {role: 'dbAdminAnyDatabase', db: 'admin'}
  ]
})

use cerebro
db.createCollection('users')
db.getCollection('users').insertOne(
  {
    id: '$UUID',
    name: '$CEREBRO_ADMIN_NAME',
    email: '$CEREBRO_ADMIN_EMAIL',
    title: null,
    positions: ['Admin'],
    password: '$CEREBRO_ADMIN_HASHED_PASSWORD',
    roles: ['Admin', 'User', 'Report', 'Data'],
    photo: null,
    verified: true,
    created: new Date().toISOString(),
    updated: new Date().toISOString()
  }
)
EOF

# Remove last command from shell history
rm /data/db/.mongodb/mongosh/mongosh_repl_history
