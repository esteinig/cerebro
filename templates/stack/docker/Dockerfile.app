FROM node:18-bookworm AS builder

{{#if dev }}
WORKDIR /app
COPY ./app/package*.json ./
RUN npm install
{{else}}
# TODO: In production mode we can download the application directory and build from revision branch
# RUN git clone --depth=1 -b {{{ revision }}} https://github.com/esteinig/cerebro && mv cerebro/app /app
WORKDIR /app
COPY ./app ./
RUN npm install
RUN npm run build
{{/if}}

FROM node:18-bookworm

ARG UID=1000
ARG GID=1000

{{#if dev}}
ENV NODE_ENV=development

# Create a user/group that matches the host, so bind mounts are writable
RUN groupadd -g ${GID} appgroup \
 && useradd  -u ${UID} -g appgroup -m appuser

WORKDIR /home/appuser/app

# Copy deps from builder; ensure ownership for the runtime user
COPY --from=builder /app/package.json ./
COPY --chown=appuser:appgroup --from=builder /app/node_modules ./node_modules

ENV PORT=8000
ENV PROTOCOL_HEADER=x-forwarded-proto
ENV HOST_HEADER=x-forwarded-host

USER appuser
{{else}}
ENV NODE_ENV=production

WORKDIR /home/node/app
COPY --from=builder /app/package.json ./
COPY --from=builder /app/node_modules ./node_modules
COPY --from=builder /app/build ./build

ENV PORT=8000
ENV PROTOCOL_HEADER=x-forwarded-proto
ENV HOST_HEADER=x-forwarded-host

USER node
{{/if}}