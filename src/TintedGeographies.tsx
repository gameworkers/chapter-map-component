import React, { useMemo, memo } from "react";
import { Geography } from "@sux/react-simple-maps";
import { Geo, Member } from "./types";

const countryNameKeys = [
  "ABBREV",
  "CONTINENT",
  "FORMAL_EN",
  "ISO_A2",
  "ISO_A3",
  "NAME",
  "NAME_LONG",
] as const;

const geographyMatchesCountryString = (
  geography: Geo,
  countryString?: string
) => countryNameKeys.some((key) => geography.properties[key] === countryString);

const TintedGeographies = ({
  geo,
  zoom,
  filter,
  memberData,
}: {
  geo: Geo[];
  zoom: number;
  filter: (geography: Geo) => boolean;
  memberData: Member[];
}) => {
  const { hasMember, noMember } = useMemo(() => {
    /* eslint-disable no-shadow */
    const hasMember: Geo[] = [];
    const noMember: Geo[] = [];
    /* eslint-enable no-shadow */

    for (const g of geo.filter(filter)) {
      if (
        memberData.some((member) =>
          geographyMatchesCountryString(g, member.country)
        )
      ) {
        hasMember.push(g);
      } else {
        noMember.push(g);
      }
    }

    return {
      hasMember,
      noMember,
    };
  }, [filter, geo, memberData]);

  const memberStyle = {
    fill: "url(#redpattern)",
    stroke: "#222",
    strokeWidth: 0.5 / zoom,
    outline: "none",
  };

  const normalStyle = {
    fill: "url(#hardlyredpattern)",
    outline: "none",
  };

  return (
    <>
      {...hasMember.map((g) => (
        <Geography
          key={g.properties["NAME"]}
          geography={g}
          style={{
            default: memberStyle,
            hover: memberStyle,
            pressed: memberStyle,
          }}
        />
      ))}
      {...noMember.map((g) => (
        <Geography
          key={g.properties["NAME"]}
          geography={g}
          style={{
            default: normalStyle,
            hover: normalStyle,
            pressed: normalStyle,
          }}
        />
      ))}
    </>
  );
};

export default memo(TintedGeographies);
