"""Collections of configured protein sequence channels."""

import hashlib
import importlib
import json

from sieve.channels import Channel, ChannelBackground, iter_fasta_records


BACKGROUND_FILE_TYPE = "sieve-channel-backgrounds"
BACKGROUND_SCHEMA_VERSION = 1
RESERVED_COLUMNS = frozenset({"sequence_id", "position", "residue"})


class ChannelSet:

    def __init__(self, channels):
        self.channels = tuple(channels)
        if not self.channels:
            raise ValueError("ChannelSet must contain at least one channel")
        for channel in self.channels:
            if not isinstance(channel, Channel):
                raise TypeError("ChannelSet entries must be Channel instances")
            if channel.short_name in RESERVED_COLUMNS:
                raise ValueError(
                    f"channel short name {channel.short_name!r} is reserved"
                )
        short_names = [channel.short_name for channel in self.channels]
        duplicates = sorted(
            {name for name in short_names if short_names.count(name) > 1}
        )
        if duplicates:
            raise ValueError(
                "duplicate channel short names: " + ", ".join(duplicates)
            )

    def __iter__(self):
        return iter(self.channels)

    def __len__(self):
        return len(self.channels)

    def definitions(self):
        return [channel.definition() for channel in self.channels]

    def fingerprint(self):
        serialized = json.dumps(
            self.definitions(),
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        return hashlib.sha256(serialized.encode("utf-8")).hexdigest()

    def compute_backgrounds(
        self,
        fasta_path,
        on_record=None,
        on_invalid_sequence=None,
    ):
        """Compute all required backgrounds in one pass over ``fasta_path``."""
        accumulators = {
            channel.short_name: channel.background_accumulator()
            for channel in self.channels
            if channel.requires_background
        }
        for sequence_id, sequence in iter_fasta_records(
            fasta_path,
            on_invalid_sequence=on_invalid_sequence,
        ):
            if on_record is not None:
                on_record(sequence_id)
            for channel in self.channels:
                accumulator = accumulators.get(channel.short_name)
                if accumulator is None:
                    continue
                if channel.background_sequence_is_usable(sequence):
                    accumulator.add_sequence(sequence)
        backgrounds = {}
        for channel in self.channels:
            accumulator = accumulators.get(channel.short_name)
            backgrounds[channel.short_name] = (
                accumulator.result() if accumulator is not None else None
            )
        return backgrounds

    def _validated_backgrounds(self, backgrounds):
        if not isinstance(backgrounds, dict):
            raise TypeError("backgrounds must be a dictionary")
        expected = {channel.short_name for channel in self.channels}
        actual = set(backgrounds)
        if actual != expected:
            missing = sorted(expected - actual)
            extra = sorted(actual - expected)
            details = []
            if missing:
                details.append("missing " + ", ".join(missing))
            if extra:
                details.append("unexpected " + ", ".join(extra))
            raise ValueError(
                "background channels do not match ChannelSet ("
                + "; ".join(details)
                + ")"
            )
        for channel in self.channels:
            background = backgrounds[channel.short_name]
            if channel.requires_background:
                if not isinstance(background, ChannelBackground):
                    raise TypeError(
                        f"channel {channel.short_name!r} requires a ChannelBackground"
                    )
            elif background is not None:
                raise ValueError(
                    f"channel {channel.short_name!r} must not have a background"
                )
        return backgrounds

    def make_functions(self, backgrounds):
        backgrounds = self._validated_backgrounds(backgrounds)
        return {
            channel.short_name: channel.make_function(
                backgrounds[channel.short_name]
            )
            for channel in self.channels
        }

    def save_backgrounds(self, backgrounds, output_path):
        backgrounds = self._validated_backgrounds(backgrounds)
        document = {
            "file_type": BACKGROUND_FILE_TYPE,
            "schema_version": BACKGROUND_SCHEMA_VERSION,
            "channel_set_fingerprint": self.fingerprint(),
            "channels": [
                {
                    "definition": channel.definition(),
                    "background": (
                        backgrounds[channel.short_name].as_dict()
                        if backgrounds[channel.short_name] is not None
                        else None
                    ),
                }
                for channel in self.channels
            ],
        }
        with open(output_path, "w", encoding="utf-8") as output:
            json.dump(document, output, indent=2, sort_keys=True, allow_nan=False)
            output.write("\n")

    def load_backgrounds(self, background_path):
        with open(background_path, encoding="utf-8") as background_file:
            document = json.load(background_file)
        if not isinstance(document, dict):
            raise ValueError("background file must contain a JSON object")
        if document.get("file_type") != BACKGROUND_FILE_TYPE:
            raise ValueError("background file has an unsupported file type")
        if document.get("schema_version") != BACKGROUND_SCHEMA_VERSION:
            raise ValueError("background file has an unsupported schema version")
        if document.get("channel_set_fingerprint") != self.fingerprint():
            raise ValueError(
                "background file does not match the configured ChannelSet"
            )
        rows = document.get("channels")
        if not isinstance(rows, list) or len(rows) != len(self.channels):
            raise ValueError("background file has invalid channel metadata")
        backgrounds = {}
        for channel, row in zip(self.channels, rows):
            if (
                not isinstance(row, dict)
                or row.get("definition") != channel.definition()
            ):
                raise ValueError(
                    "background file channel definitions do not match ChannelSet"
                )
            raw_background = row.get("background")
            if channel.requires_background:
                if not isinstance(raw_background, dict):
                    raise ValueError(
                        f"background for channel {channel.short_name!r} is missing"
                    )
                try:
                    background = ChannelBackground(
                        raw_background["bg_mu"], raw_background["bg_sigma"]
                    )
                except KeyError as error:
                    raise ValueError(
                        f"background for channel {channel.short_name!r} is incomplete"
                    ) from error
            else:
                if raw_background is not None:
                    raise ValueError(
                        f"channel {channel.short_name!r} must not have a background"
                    )
                background = None
            backgrounds[channel.short_name] = background
        return self._validated_backgrounds(backgrounds)


def load_channel_set(spec):
    """Load a ChannelSet from a dotted ``module.attribute`` path."""
    module_name, separator, attribute_name = spec.rpartition(".")
    if not separator:
        raise ValueError("ChannelSet spec must be in the form module.attribute")
    module = importlib.import_module(module_name)
    try:
        channel_set = getattr(module, attribute_name)
    except AttributeError as error:
        raise ValueError(f"Cannot find ChannelSet attribute: {spec}") from error
    if not isinstance(channel_set, ChannelSet):
        raise TypeError(
            f"{spec} must refer to a sieve.channel_set.ChannelSet instance"
        )
    return channel_set
