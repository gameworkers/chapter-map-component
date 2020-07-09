import useSWR from "swr";

/**
 * swr-based fetcher configured to fetch only once. still helps handle some
 * control flow issues like component unmounting, etc., over just using
 * useEffect + fetch + useState.
 */
const useFetch = <T>(key: string) =>
  useSWR<T>(key, (url: string) => fetch(url).then((res) => res.json()), {
    revalidateOnFocus: false,
    revalidateOnReconnect: false,
    refreshWhenOffline: false,
    refreshWhenHidden: false,
    refreshInterval: 0,
  });

export default useFetch;
