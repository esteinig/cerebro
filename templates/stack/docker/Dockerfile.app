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

ENV PORT=8000
ENV PROTOCOL_HEADER=x-forwarded-proto
ENV HOST_HEADER=x-forwarded-host

# single path for both dev/prod
WORKDIR /app

{{#if dev}}
ENV NODE_ENV=development

# Create/reuse UID/GID safely (handles UID/GID=1000 already occupied by node)
RUN set -eux; \
    if ! getent group "${GID}" >/dev/null; then groupadd -g "${GID}" appgroup; fi; \
    if ! getent passwd "${UID}" >/dev/null; then useradd -u "${UID}" -g "${GID}" -m appuser; fi

# copy deps with numeric ownership so we don't care about usernames
COPY --from=builder /app/package.json ./
COPY --from=builder /app/node_modules ./node_modules
RUN chown -R ${UID}:${GID} /app

USER ${UID}:${GID}
{{else}}
ENV NODE_ENV=production

COPY --from=builder /app/package.json ./
COPY --from=builder /app/node_modules ./node_modules
COPY --from=builder /app/build ./build

USER node
{{/if}}