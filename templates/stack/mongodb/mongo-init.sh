set -e

UUID=$(cat /proc/sys/kernel/random/uuid)

# Service account identifiers (S3-5 #5): a dedicated Bot user and a service Team
# that owns its lifecycle data. The Faktory lifecycle worker authenticates as this
# Bot (role-scoped, long-lived token) and operates within this team's admin
# database. All identifiers are generated here so the team's `users` array and the
# admin database name line up across the documents.
BOT_UUID=$(cat /proc/sys/kernel/random/uuid)
TEAM_UUID=$(cat /proc/sys/kernel/random/uuid)
ADMIN_DB_UUID=$(cat /proc/sys/kernel/random/uuid)
USER_DB_UUID=$(cat /proc/sys/kernel/random/uuid)
USER_PROJECT_UUID=$(cat /proc/sys/kernel/random/uuid)
P_FILES=$(cat /proc/sys/kernel/random/uuid)
P_WATCHERS=$(cat /proc/sys/kernel/random/uuid)
P_PIPELINES=$(cat /proc/sys/kernel/random/uuid)
P_LOGS=$(cat /proc/sys/kernel/random/uuid)
P_REPORTS=$(cat /proc/sys/kernel/random/uuid)
P_TRAINDATA=$(cat /proc/sys/kernel/random/uuid)
P_TRAINSESS=$(cat /proc/sys/kernel/random/uuid)

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

// Service Bot user for the Faktory lifecycle worker (S3-5 #5). Must hold the
// 'Bot' role (the API only mints a long-lived bot token via ?role=Bot when the
// user has it) and 'User'/'Data' for team-scoped data access. Verified so it can
// authenticate immediately without the email verification flow.
db.getCollection('users').insertOne(
  {
    id: '$BOT_UUID',
    name: '$CEREBRO_SERVICE_BOT_NAME',
    email: '$CEREBRO_SERVICE_BOT_EMAIL',
    title: null,
    positions: ['Service'],
    password: '$CEREBRO_SERVICE_BOT_HASHED_PASSWORD',
    roles: ['Bot', 'User', 'Data'],
    photo: null,
    verified: true,
    created: new Date().toISOString(),
    updated: new Date().toISOString()
  }
)

// Service Team owned by the Bot. The worker addresses files via ?team=<this team>;
// they live in the team's admin database (admin.database below). The admin
// database's project set mirrors what the API would create for a team, so this
// team is indistinguishable from one created through the application.
db.createCollection('teams')
db.getCollection('teams').insertOne(
  {
    id: '$TEAM_UUID',
    name: '$CEREBRO_SERVICE_TEAM_NAME',
    description: 'Service team for the Cerebro lifecycle worker',
    users: ['$BOT_UUID'],
    databases: [
      {
        id: '$USER_DB_UUID',
        name: '$CEREBRO_SERVICE_DB_NAME',
        database: '$USER_DB_UUID',
        description: 'Service database',
        projects: [
          {
            id: '$USER_PROJECT_UUID',
            name: '$CEREBRO_SERVICE_PROJECT_NAME',
            collection: '$USER_PROJECT_UUID',
            description: 'Default data collection'
          }
        ]
      }
    ],
    admin: {
      id: '$ADMIN_DB_UUID',
      name: 'Admin',
      database: '$ADMIN_DB_UUID',
      description: 'Team administrative database',
      projects: [
        { id: '$P_FILES',     name: 'Files',             collection: 'files',             description: 'File registrations for CerebroFS' },
        { id: '$P_WATCHERS',  name: 'Watchers',          collection: 'watchers',          description: 'Production watcher registrations for file uploads to CerebroFS' },
        { id: '$P_PIPELINES', name: 'Pipelines',         collection: 'towers',            description: 'Production pipeline registrations for Cerebro' },
        { id: '$P_LOGS',      name: 'Logs',              collection: 'logs',              description: 'Log entries for the team' },
        { id: '$P_REPORTS',   name: 'Reports',           collection: 'reports',           description: 'Report storage for team independent of database' },
        { id: '$P_TRAINDATA', name: 'TrainingData',      collection: 'training_data',     description: 'Training data for the team' },
        { id: '$P_TRAINSESS', name: 'TrainingSessions',  collection: 'training_sessions', description: 'Training sessions for the team' }
      ]
    }
  }
)
EOF

# Remove last command from shell history
rm /data/db/.mongodb/mongosh/mongosh_repl_history
