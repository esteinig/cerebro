import type { PageServerLoad } from './$types';
import { postJSON } from '$lib/server/api';
import { env as private_env } from "$env/dynamic/private";

type VerificationCheck =
  | { status: 'ok'; message?: string; access_token: string }
  | { status: 'already_used'; message?: string }

export const load: PageServerLoad = async ({ url }) => {
  const token = url.searchParams.get('token');
  if (!token) return { step: 'needs_token' as const };

  try {
    const check = await postJSON<VerificationCheck>(
      '/auth/verification-check',
      { access_token: token },
      private_env.PRIVATE_CEREBRO_API_URL_DOCKER
    );

    if (check.status === 'ok') {
      await postJSON(
        '/auth/password-reset-email',
        { access_token: check.access_token },
        private_env.PRIVATE_CEREBRO_API_URL_DOCKER
      );
      return { step: 'email_sent' as const };
    }

    if (check.status === 'already_used') {
      // Idempotent success: treat as sent
      return { step: 'email_sent' as const };
    }

    return { step: 'invalid' as const };
  } catch {
    return { step: 'invalid' as const };
  }
};