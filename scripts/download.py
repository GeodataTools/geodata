import argparse
import logging
import os

logging.basicConfig(level=logging.INFO)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--save-dir")
    parser.add_argument("--year-start", "-ys", required=True, type=int)
    parser.add_argument("--year-end", "-ye", required=True, type=int)
    parser.add_argument("--month-start", "-ms", type=int)
    parser.add_argument("--month-end", "-me", type=int)
    parser.add_argument("--config", "-c", required=True)
    parser.add_argument("--module", "-m", required=True)
    args = parser.parse_args()

    if args.month_start is None or args.month_end is None:
        months = slice(1, 12)
    else:
        months = slice(args.month_start, args.month_end)

    if args.save_dir is not None:
        os.environ["GEODATA_ROOT"] = args.save_dir
    import geodata

    dataset = geodata.Dataset(
        module=args.module,
        weather_data_config=args.config,
        years=slice(args.year_start, args.year_end),
        months=months,
    )

    if not dataset.prepared:
        dataset.get_data()


if __name__ == "__main__":
    main()
