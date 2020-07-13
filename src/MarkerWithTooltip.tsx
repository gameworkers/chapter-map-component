import React, { memo } from "react";

import TooltipTrigger, { TooltipArg, ChildrenArg } from "react-popper-tooltip";

import GWUMarker from "./GWUMarker";

import type { Marker } from "./types";

const MarkerWithTooltip = ({
  marker,
  scale,
  className,
  onTooltipChange,
  showTooltip,
}: {
  marker: Marker;
  scale: number;
  className: string;
  onTooltipChange: (visible: boolean) => void;
  showTooltip: boolean;
}) => (
  <TooltipTrigger
    placement="top"
    trigger="click"
    onVisibilityChange={onTooltipChange}
    tooltipShown={showTooltip}
    closeOnReferenceHidden
    tooltip={(args: TooltipArg) => (
      <Tooltip {...args}>
        <MarkerTooltip marker={marker} className={className} />
      </Tooltip>
    )}
  >
    {({ triggerRef, getTriggerProps }: ChildrenArg) => (
      <GWUMarker
        ref={triggerRef as any}
        marker={marker}
        scale={scale}
        {...getTriggerProps()}
      />
    )}
  </TooltipTrigger>
);

const Tooltip = ({
  arrowRef,
  tooltipRef,
  getArrowProps,
  getTooltipProps,
  placement,
  children,
}: TooltipArg & {
  children: React.ReactNode;
}) => (
  <div
    {...getTooltipProps({
      ref: tooltipRef,
      className: "tooltip-container",
    })}
  >
    <div
      {...getArrowProps({
        ref: arrowRef,
        className: "tooltip-arrow",
        "data-placement": placement,
      })}
    />
    {children}
  </div>
);

const MarkerTooltip = ({
  marker: {
    data: { location, chapterInfo },
  },
  className,
}: {
  marker: Marker;
  className?: string;
}) => {
  if (!chapterInfo) return null;
  const { description, applicationLink, twitter, email, website } = chapterInfo;
  return (
    <div className={className}>
      <h3>{location}</h3>
      {description && <p dangerouslySetInnerHTML={{ __html: description }} />}
      {(twitter || email || website) && (
        <p>
          {twitter && (
            <>
              <strong>Twitter: </strong>
              <a href={`https://twitter.com/${twitter}`}>@{twitter}</a>
              <br />
            </>
          )}
          {email && (
            <>
              <strong>Email: </strong>
              <a href={`mailto:${email}`}>{email}</a>
              <br />
            </>
          )}
          {website && (
            <>
              <strong>Website: </strong>
              <a
                href={`${
                  website.startsWith("http") ? "" : "http://"
                }${website}`}
              >
                {website}
              </a>
              <br />
            </>
          )}
        </p>
      )}
      {applicationLink && (
        <p>
          <a href={applicationLink}>Apply here!</a>
        </p>
      )}
    </div>
  );
};

export default memo(MarkerWithTooltip);
