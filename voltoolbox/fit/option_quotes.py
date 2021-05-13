import json
import dataclasses
from dataclasses import dataclass
from typing import Union, Dict, List, Tuple
from datetime import date, datetime

class EnhancedJSONEncoder(json.JSONEncoder):
    def default(self, o):
        if dataclasses.is_dataclass(o):
            return dataclasses.asdict(o)
        elif isinstance(o, (date, datetime)):
            return o.isoformat()
        else:
            return super().default(o)


class JsonableObject:
    
    def to_json(self) -> str:
        return json.dumps(self, cls=EnhancedJSONEncoder, indent=4)


@dataclass(frozen=True)
class QuoteSlice(JsonableObject):
    strikes: Tuple[float]
    bids: Tuple[float]
    asks: Tuple[float]

    @classmethod
    def from_json_dict(cls, json_dict):
        return cls(tuple(json_dict['strikes']),
                   tuple(json_dict['bids']),
                   tuple(json_dict['asks']))


@dataclass(frozen=True)
class OptionQuoteSlice(JsonableObject):
    symbol: str
    expiry: date    
    call : QuoteSlice
    put : QuoteSlice

    @classmethod
    def from_json_dict(cls, json_dict):
        return cls(str(json_dict['symbol']),
                   date.fromisoformat(json_dict['expiry']),
                   QuoteSlice.from_json_dict(json_dict['call']),
                   QuoteSlice.from_json_dict(json_dict['put']))


@dataclass(frozen=True)
class OptionSnapshot(JsonableObject):
    time_stamp: datetime
    ref_spot: float
    slices: List[OptionQuoteSlice]

    @classmethod
    def from_json_dict(cls, json_dict):
        slices = [OptionQuoteSlice.from_json_dict(sl) 
                  for sl in json_dict['slices']]
        return cls(datetime.fromisoformat(json_dict['time_stamp']),
                   float(json_dict['ref_spot']),
                   slices)

