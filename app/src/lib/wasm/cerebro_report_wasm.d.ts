/* tslint:disable */
/* eslint-disable */
export function slugify(s: string): string;
export class ReportCompiler {
  free(): void;
  constructor(root: string, request_data: Function, report_type: string, report_logo?: Uint8Array);
  report(config: any): string;
  pdf(text: string, path: string): Uint8Array;
  svg(text: string, path: string): (string)[];
  add_font(data: Uint8Array): void;
  add_logo(data: Uint8Array): void;
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
  readonly memory: WebAssembly.Memory;
  readonly __wbg_reportcompiler_free: (a: number, b: number) => void;
  readonly reportcompiler_new: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => void;
  readonly reportcompiler_report: (a: number, b: number, c: number) => void;
  readonly reportcompiler_pdf: (a: number, b: number, c: number, d: number, e: number, f: number) => void;
  readonly reportcompiler_svg: (a: number, b: number, c: number, d: number, e: number, f: number) => void;
  readonly reportcompiler_add_font: (a: number, b: number, c: number) => void;
  readonly reportcompiler_add_logo: (a: number, b: number, c: number, d: number) => void;
  readonly slugify: (a: number, b: number, c: number) => void;
  readonly qcms_transform_data_rgb_out_lut_precache: (a: number, b: number, c: number, d: number) => void;
  readonly qcms_transform_data_rgba_out_lut_precache: (a: number, b: number, c: number, d: number) => void;
  readonly qcms_transform_data_bgra_out_lut_precache: (a: number, b: number, c: number, d: number) => void;
  readonly qcms_transform_data_rgb_out_lut: (a: number, b: number, c: number, d: number) => void;
  readonly qcms_transform_data_rgba_out_lut: (a: number, b: number, c: number, d: number) => void;
  readonly qcms_transform_data_bgra_out_lut: (a: number, b: number, c: number, d: number) => void;
  readonly qcms_transform_release: (a: number) => void;
  readonly qcms_profile_precache_output_transform: (a: number) => void;
  readonly qcms_enable_iccv4: () => void;
  readonly qcms_profile_is_bogus: (a: number) => number;
  readonly qcms_white_point_sRGB: (a: number) => void;
  readonly lut_interp_linear16: (a: number, b: number, c: number) => number;
  readonly lut_inverse_interp16: (a: number, b: number, c: number) => number;
  readonly __wbindgen_malloc: (a: number, b: number) => number;
  readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
  readonly __wbindgen_exn_store: (a: number) => void;
  readonly __wbindgen_add_to_stack_pointer: (a: number) => number;
  readonly __wbindgen_free: (a: number, b: number, c: number) => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;
/**
* Instantiates the given `module`, which can either be bytes or
* a precompiled `WebAssembly.Module`.
*
* @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
*
* @returns {InitOutput}
*/
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
* If `module_or_path` is {RequestInfo} or {URL}, makes a request and
* for everything else, calls `WebAssembly.instantiate` directly.
*
* @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
*
* @returns {Promise<InitOutput>}
*/
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
