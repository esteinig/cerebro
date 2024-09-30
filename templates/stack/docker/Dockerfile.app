FROM node:18-bookworm AS builder

{{#if dev }}
# In development mode we mount the application path into the container
WORKDIR /app
COPY ./app/package*.json ./
RUN npm install
{{else}}
# In production mode we download the application directory and build from revision branch
# RUN git clone --depth=1 -b {{{ revision }}} https://github.com/esteinig/cerebro && mv cerebro/app /app
WORKDIR /app
COPY ./app ./
RUN npm install
RUN npm run build
{{/if}}

FROM node:18-bookworm

{{#if dev }}
ENV NODE_ENV=development
{{else}}
ENV NODE_ENV=production
{{/if}}

WORKDIR /home/node/app
COPY --from=builder /app/package.json .

{{#if dev}}
COPY --chown=node:node --from=builder /app/node_modules ./node_modules
{{else}}
COPY --from=builder /app/node_modules ./node_modules
COPY --from=builder /app/build ./build
{{/if}}

ENV PORT=8000
ENV PROTOCOL_HEADER=x-forwarded-proto
ENV HOST_HEADER=x-forwarded-host 

USER node
