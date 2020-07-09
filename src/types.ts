/**
 * Models data retrieved from member data JSON file
 */
export interface Member {
  isUnion?: boolean;
  isChapter?: boolean;
  lat: number;
  lng: number;
  location: string;
  country?: string;
  chapterInfo: {
    description?: string;
    applicationLink?: string;
    twitter?: string;
    email?: string;
    website?: string;
  };
}

export interface Marker {
  name: string;
  coordinates: [number, number];
  data: Member;
}

/**
 * Models data retrieved from geodata JSON, eg world-50m.json
 */
export interface Geo {
  properties: {
    [key: string]: string;
  };
}
