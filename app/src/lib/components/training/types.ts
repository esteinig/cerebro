import { AbundanceMode, DisplayData, DisplayVisualisation} from "$lib/utils/types";
import type { TaxonFilterConfig } from "$lib/utils/types";

export type UUID = string;

export enum SlideDecision {
  Positive = "Positive",
  Negative = "Negative"
}

export type SlideId = string;

export type SlideConfig = {
  id: SlideId;
  title: string;
  description?: string;
  displayData: DisplayData;
  displayMode: AbundanceMode;
  selectedVisualisation: DisplayVisualisation; // typically DisplayVisualisation.Table or Threshold
  serverFilterConfig?: TaxonFilterConfig | null;
  disablePrevalenceContamination?: boolean;
  // Optional hard overrides
  forceRpmForValidationPlate?: boolean;
}

export type DeckConfig = {
  id: UUID;
  name: string;
  slides: SlideConfig[];
  minSelection?: number;        // default 1
  maxSelection?: number;        // optional
}

export type SlideResponse = {
  slideId: SlideId;
  decision: SlideDecision;
  selectedTaxa: { taxid: string; name: string }[]; // snapshot of chosen rows
  msSpent: number;                                  // time spent on this slide
  comment?: string;
  clientEventId: string;                            // to prevent duplicates
}

export type TrainingSessionStart = {
  deckId: UUID;
  userId: UUID;
  sessionId: UUID;
  startedAt: string; // ISO
}

export type TrainingSessionComplete = {
  sessionId: UUID;
  totalMs: number;
  slideCount: number;
  responses: SlideResponse[]; // echoed back or returned on GET
}
