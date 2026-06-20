

from og.core.generic import Enumerative


ENUM_CALLKEY = Enumerative(
    INCLUDE_GROUPS='include-groups',
    EXCLUDE_GROUPS='exclude-groups',
    INCLUDE_MEMBERS='include-members',
    EXCLUDE_MEMBERS='exclude-members',
    INCLUDE_SAMPLES='include-samples',
    EXCLUDE_SAMPLES='exclude-samples'
)


def subset(ortho, call_order):
    call_function = dict()
    def call_include_groups(include_groups):
        groups = []
        for group in ortho.groups:
            if group.id in include_groups:
                groups.append(group)
        ortho.groups = groups
    call_function[ENUM_CALLKEY.INCLUDE_GROUPS] = call_include_groups

    def call_exclude_groups(exclude_groups):
        groups = []
        for group in ortho.groups:
            if group.id not in exclude_groups:
                groups.append(group)
        ortho.groups = groups
    call_function[ENUM_CALLKEY.EXCLUDE_GROUPS] = call_exclude_groups

    def call_include_samples(include_samples):
        ortho.subset_samples(include_samples, inplace=True)
    call_function[ENUM_CALLKEY.INCLUDE_SAMPLES] = call_include_samples

    def call_exclude_samples(exclude_samples):
        include_samples = []
        for sample in ortho.samples:
            if sample.id not in exclude_samples:
                include_samples.append(sample.id)
        ortho.subset_samples(include_samples, inplace=True)
    call_function[ENUM_CALLKEY.EXCLUDE_SAMPLES] = call_exclude_samples

    def call_include_members(include_members):
        for group in ortho.groups:
            for sample in ortho.samples:
                members = []
                for member in group[sample.index]:
                    if member in include_members:
                        members.append(member)
                group[sample.index] = members
    call_function[ENUM_CALLKEY.INCLUDE_MEMBERS] = call_include_members

    def call_exclude_members(exclude_members):
        for group in ortho.groups:
            for sample in ortho.samples:
                members = []
                for member in group[sample.index]:
                    if member not in exclude_members:
                        members.append(member)
                group[sample.index] = members
    call_function[ENUM_CALLKEY.EXCLUDE_MEMBERS] = call_exclude_members
    
    for call_key, call_input in call_order:
        call_function[call_key](call_input)

    return ortho
