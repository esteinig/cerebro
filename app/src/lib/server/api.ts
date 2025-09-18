// Server-only helper for calling API
export async function postJSON<T>(
  path: string,
  body: unknown,
  api_base: string,
  signal?: AbortSignal
): Promise<T> {
  const res = await fetch(`${api_base}${path}`, {
    method: 'POST',
    headers: { 'content-type': 'application/json' },
    body: JSON.stringify(body),
    signal
  });
  if (!res.ok) {
    throw new Error(`Request failed with  status: ${res.status}`);
  }
  return (await res.json()) as T;
}


export function parseBool(value: string | undefined): boolean {
    if (value === 'true') return true;
    if (value === 'false') return false;
    throw new Error(`Invalid Configuration`);
}

export function parseSameSite(value: string | undefined): 'lax' | 'strict' | 'none' | undefined {
    if (!value) return undefined;
    const lower = value.toLowerCase();
    if (lower === 'lax' || lower === 'strict' || lower === 'none') return lower;
    throw new Error(`Invalid SameSite value: ${value}`);
  }
  