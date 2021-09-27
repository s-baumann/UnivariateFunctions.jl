const tol = 10*eps()
const days_per_year = 365.2422
const global_base_date = Date(2000,1,1)
const global_base_date_as_day = convert(Dates.Day, global_base_date)

"""
    years_between(a::Union{DateTime,Date}, b::Union{DateTime,Date})
The number of years between two dates. This is returned as a scalar and
assumes 365.2422 days per year.

### Inputs
* `a` - The end date
* `b` - The start date.
### Returns
* A `Real`.
"""
function years_between(a::Union{DateTime,Date}, b::Union{DateTime,Date})
    return (Dates.days(a) -Dates.days(b))/ days_per_year
end

"""
    years_from_global_base(a::Union{DateTime,Date})
The number of years (calculated by the `years_between` function) between a given date
and the 1st of January 2000.

### Inputs
* `date` - A date.
### Returns
* A `Real`.
"""
function years_from_global_base(date::Union{DateTime,Date})
    return years_between(date, global_base_date)
end


"""
    period_length(a::Dates.DatePeriod, base::Date = global_base_date)

Period length is designed to convert `TimePeriod` objects to a float in a consistent way to `years_from_global_base`.
So effectively the years_between method is calculated with start and end dates being those
at the start and end of a `Dates.DatePeriod`. This is slightly complicated because
a period like `Month(3)` might have slightly different numbers of total days depending
on when in the year it is. So a base date has to be input. The period is then
measured starting from this base date.

### Inputs
* `period` - A period.
* `base` - A date from which the period will be measured from.
### Returns
* A `Real`.
"""
function period_length(period::Dates.DatePeriod, base::Date = global_base_date)
    return years_between(base+period, base)
end
