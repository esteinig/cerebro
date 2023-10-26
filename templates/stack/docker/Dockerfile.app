{{#if dev }}
ENV NODE_ENV=development
{{else}}
ENV NODE_ENV=production
{{/if}}

FROM node:18-bookworm AS builder

{{#if dev }}
WORKDIR /app
COPY ./app/package*.json ./
RUN npm install
# In development mode we mount the application path into the container
{{else}}
WORKDIR /src
# In production mode we download the application directory and build from revision branch
RUN git clone --depth=1 -b {{{ revision }}} https://github.com/esteinig/cerebro && mv cerebro/app /app
WORKDIR /app
RUN npm install
RUN npm run build
{{/if}}


FROM node:18-bookworm

WORKDIR /home/node/app
COPY --from=builder /app/package.json .
COPY --from=builder /app/node_modules ./node_modules
{{#unless dev }}
COPY --from=builder /app/build ./build
{{/unless}}

ENV PORT=8000
ENV PROTOCOL_HEADER=x-forwarded-proto
ENV HOST_HEADER=x-forwarded-host 

USER node
